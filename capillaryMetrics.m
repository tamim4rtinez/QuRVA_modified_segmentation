function [stats, qc] = capillaryMetrics(im_clean, im_skel, params)
% capillaryMetrics  Morphometric analysis of the retinal capillary bed.
%
%   stats        = capillaryMetrics()                      % NaN template
%   [stats, qc]  = capillaryMetrics(im_clean, im_skel, params)
%
% Separates capillaries from large vessels by local vessel radius, then
% quantifies complexity, branch geometry, avascular spaces, perfusion
% distance, caliber, and spatial heterogeneity.
%
% params fields (all optional except the first two):
%   radiusThreshPix  local RADIUS cutoff, = HALF the vessel width. Default 6.
%   minBranchPix     discard branches shorter than this as spurs. Default 8.
%   pixSzeUM         microns per pixel. Default 1 (outputs then in pixels).
%   hull             leaflet-area mask. Default bwconvhull(im_clean).
%   gridN            tiles per side for spatial CV. Default 8.
%   apexXY           [x y] of the optic-nerve apex, for radial binning.
%                    If omitted, radial metrics return NaN.
%   nRadialBins      number of radial bins. Default 6.

    % ------------------------------------------------------------------
    % Template mode: return the scalar field list with NaNs
    % ------------------------------------------------------------------
    if nargin == 0
        stats = emptyStats();
        return
    end

    if nargin < 3, params = struct(); end
    p = defaultParams(params);

    im_clean = logical(im_clean);
    im_skel  = logical(im_skel);

    if isempty(p.hull)
        im_hull = bwconvhull(im_clean, 'union');
    else
        im_hull = logical(p.hull);
    end

    px  = p.pixSzeUM;
    px2 = px^2;

    stats = emptyStats();

    % ==================================================================
    % 1. Split skeleton by local vessel radius
    % ==================================================================
    D = bwdist(~im_clean);          % at a skeleton point this IS the radius

    skelRadius = zeros(size(D));
    skelRadius(im_skel) = D(im_skel);

    skel_cap   = im_skel & (skelRadius > 0) & (skelRadius <= p.radiusThreshPix);
    skel_large = im_skel & (skelRadius > p.radiusThreshPix);

    skel_cap = bwareaopen(skel_cap, p.minBranchPix, 8);

    % ==================================================================
    % 2. Junctions, endpoints, components
    % ==================================================================
    nbrCount    = conv2(double(skel_cap), [1 1 1; 1 0 1; 1 1 1], 'same');
    junctionPix = skel_cap & (nbrCount >= 3);
    endPix      = skel_cap & (nbrCount == 1);

    ccJunc     = bwconncomp(junctionPix, 8);
    nJunctions = ccJunc.NumObjects;
    nEndpoints = nnz(endPix);
    nComponents = bwconncomp(skel_cap, 8).NumObjects;

    if nJunctions > 0
        jc = regionprops(ccJunc, 'Centroid');
        juncXY = cat(1, jc.Centroid);       % [x y]
    else
        juncXY = zeros(0,2);
    end

    % ==================================================================
    % 3. Branch lengths (true path length) and per-branch tortuosity
    % ==================================================================
    branches = skel_cap & ~junctionPix;
    L        = bwlabel(branches, 8);
    nLabels  = max(L(:));

    lengths = zeros(nLabels, 1);
    offsets = {[0 1], [1 0], [1 1], [1 -1]};
    stepLen = [1, 1, sqrt(2), sqrt(2)];

    for k = 1:4
        Ls = circshift(L, offsets{k});
        if offsets{k}(1) ~= 0,  Ls(1,:)   = 0; end
        if offsets{k}(2) ==  1, Ls(:,1)   = 0; end
        if offsets{k}(2) == -1, Ls(:,end) = 0; end
        same = (L > 0) & (L == Ls);
        if any(same(:))
            lengths = lengths + accumarray(L(same), stepLen(k), [nLabels 1]);
        end
    end

    pixCounts = accumarray(L(L>0), 1, [nLabels 1]);
    lengths(lengths == 0 & pixCounts > 0) = 1;

    keep       = lengths >= p.minBranchPix;
    lengths_um = lengths(keep) * px;

    % --- tortuosity: path length / straight-line distance between the two
    %     free ends of each branch. Branches forming closed loops (0 or >2
    %     endpoints) are skipped rather than given a meaningless ratio.
    brNbr  = conv2(double(branches), [1 1 1; 1 0 1; 1 1 1], 'same');
    epMask = branches & (brNbr <= 1);
    epLab  = L(epMask);
    [epR, epC] = find(epMask);

    tort = nan(nLabels, 1);
    if ~isempty(epLab)
        epCount = accumarray(epLab, 1, [nLabels 1]);
        twoEnds = find(epCount == 2 & lengths > 0);
        for m = 1:numel(twoEnds)
            lb  = twoEnds(m);
            idx = find(epLab == lb);
            dEuc = hypot(epR(idx(1)) - epR(idx(2)), epC(idx(1)) - epC(idx(2)));
            if dEuc > 0
                tort(lb) = lengths(lb) / dEuc;
            end
        end
    end
    tort = tort(keep);
    tort = tort(isfinite(tort) & tort >= 1);

    % ==================================================================
    % 4. Loops in the capillary network (web vs tree)
    % ==================================================================
    holeCount     = nComponents - bweuler(skel_cap, 8);
    holeCountMask = bwconncomp(im_clean,8).NumObjects - bweuler(im_clean, 8);

    % ==================================================================
    % 5. Avascular space size distribution (holes in the vessel mask)
    % ==================================================================
    holesMask = imfill(im_clean, 'holes') & ~im_clean;
    ccHoles   = bwconncomp(holesMask, 8);
    if ccHoles.NumObjects > 0
        ha = regionprops(ccHoles, 'Area');
        holeAreas_um2 = [ha.Area]' * px2;
    else
        holeAreas_um2 = zeros(0,1);
    end

    hullPix         = nnz(im_hull);
    leafletArea_mm2 = hullPix * px2 / 1e6;

    % ==================================================================
    % 6. Intercapillary distance (tissue distance to nearest vessel)
    % ==================================================================
    Dout    = bwdist(im_clean);
    inTissue = im_hull & ~im_clean;
    icd_um  = Dout(inTissue) * px;

    % ==================================================================
    % 7. Capillary caliber
    % ==================================================================
    capRadius_um = D(skel_cap) * px;

    % ==================================================================
    % 8. Spatial heterogeneity
    % ==================================================================
       [gridCV_density, gridCV_junction, nTilesUsed] = ...
        gridHeterogeneity(skel_cap, im_hull, juncXY, p.gridN);

        [radialCV_density, radialProfile] = ...
        radialHeterogeneity(skel_cap, im_hull, p.apexXY, p.nRadialBins);

    % ==================================================================
    % 9. Assemble
    % ==================================================================
    capSkelLen_mm = sum(lengths(keep)) * px / 1e3;

    stats.nJunctions             = nJunctions;
    stats.nEndpoints             = nEndpoints;
    stats.nBranches              = numel(lengths_um);
    stats.nComponents            = nComponents;
    stats.leafletArea_mm2        = leafletArea_mm2;
    stats.capSkelLen_mm          = capSkelLen_mm;

    stats.junctionsPerArea_mm2   = nJunctions / max(leafletArea_mm2, eps);
    stats.junctionsPerLength_mm  = nJunctions / max(capSkelLen_mm,   eps);
    stats.endpointJunctionRatio  = nEndpoints / max(nJunctions, 1);

    stats.holeCount              = holeCount;
    stats.holesPerArea_mm2       = holeCount / max(leafletArea_mm2, eps);
    stats.holeCountMask          = holeCountMask;

    stats.meanBranchLen_um       = safeStat(@mean,   lengths_um);
    stats.sdBranchLen_um         = safeStat(@std,    lengths_um);
    stats.medianBranchLen_um     = safeStat(@median, lengths_um);
    stats.iqrBranchLen_um        = safeStat(@iqr,    lengths_um);
    stats.cvBranchLen            = stats.sdBranchLen_um / max(stats.meanBranchLen_um, eps);

    stats.meanTortuosity         = safeStat(@mean,   tort);
    stats.sdTortuosity           = safeStat(@std,    tort);
    stats.medianTortuosity       = safeStat(@median, tort);

    stats.nHoles                 = numel(holeAreas_um2);
    stats.meanHoleArea_um2       = safeStat(@mean,   holeAreas_um2);
    stats.medianHoleArea_um2     = safeStat(@median, holeAreas_um2);
    stats.sdHoleArea_um2         = safeStat(@std,    holeAreas_um2);
    stats.maxHoleArea_um2        = safeStat(@max,    holeAreas_um2);
    stats.cvHoleArea             = stats.sdHoleArea_um2 / max(stats.meanHoleArea_um2, eps);
    stats.avascularFrac          = sum(holeAreas_um2) / max(hullPix * px2, eps);

    stats.meanIntercapDist_um    = safeStat(@mean, icd_um);
    stats.p95IntercapDist_um     = safeStat(@(x) prctile(x,95), icd_um);
    stats.maxIntercapDist_um     = safeStat(@max,  icd_um);

    stats.meanCapRadius_um       = safeStat(@mean, capRadius_um);
    stats.sdCapRadius_um         = safeStat(@std,  capRadius_um);

    stats.gridCV_density         = gridCV_density;
    stats.gridCV_junctionDensity = gridCV_junction;
    stats.nGridTiles             = nTilesUsed;
    stats.radialCV_density       = radialCV_density;
        for b = 1:min(numel(radialProfile), 6)
        stats.(sprintf('radialDens_b%d', b)) = radialProfile(b);
        end

    stats.largeSkelPixels        = nnz(skel_large);
    stats.capSkelPixels          = nnz(skel_cap);

    % Non-scalar extras (not written to the CSV; use for plotting)
    stats.branchLengths_um       = lengths_um;
    stats.branchTortuosity       = tort;
    stats.holeAreas_um2          = holeAreas_um2;
    stats.radialProfile          = radialProfile;

    % ==================================================================
    % 10. QC image
    % ==================================================================
    if nargout > 1
        basei = 0.25 * double(im_clean);
        qc    = repmat(basei, [1 1 3]);
        capD  = double(bwmorph(skel_cap,   'dilate', 1));
        lrgD  = double(bwmorph(skel_large, 'dilate', 1));
        qc(:,:,2) = max(qc(:,:,2), capD);
        qc(:,:,1) = max(qc(:,:,1), lrgD);
        qc(:,:,3) = max(qc(:,:,3), lrgD);
    end
end

% ======================================================================
function p = defaultParams(params)
    p.radiusThreshPix = getf(params, 'radiusThreshPix', 6);
    p.minBranchPix    = getf(params, 'minBranchPix',    8);
    p.pixSzeUM        = getf(params, 'pixSzeUM',        1);
    p.hull            = getf(params, 'hull',            []);
    p.gridN           = getf(params, 'gridN',           8);
    p.apexXY          = getf(params, 'apexXY',          []);
    p.nRadialBins     = getf(params, 'nRadialBins',     6);
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function v = safeStat(fn, x)
    if isempty(x), v = NaN; else, v = fn(double(x)); end
end

% ----------------------------------------------------------------------
function [cvDens, cvJunc, nUsed] = gridHeterogeneity(skel_cap, im_hull, juncXY, gridN)
% Tile the hull bounding box; vessel density and junction density per tile;
% CV across tiles that are mostly inside the tissue.

    cvDens = NaN; cvJunc = NaN; nUsed = 0;

    [r, c] = find(im_hull);
    if isempty(r), return; end
    r0 = min(r); r1 = max(r); c0 = min(c); c1 = max(c);

    rEdges = round(linspace(r0, r1+1, gridN+1));
    cEdges = round(linspace(c0, c1+1, gridN+1));

    dens = []; jdens = [];

    for i = 1:gridN
        for j = 1:gridN
            rr = rEdges(i):rEdges(i+1)-1;
            cc = cEdges(j):cEdges(j+1)-1;
            if isempty(rr) || isempty(cc), continue; end

            hullTile = im_hull(rr, cc);
            tileArea = numel(hullTile);
            hullPix  = nnz(hullTile);

            % Skip tiles mostly outside the tissue - they would otherwise
            % register as "low density" and inflate the CV artificially.
            if hullPix < 0.5 * tileArea, continue; end

            skelPix = nnz(skel_cap(rr, cc) & hullTile);
            dens(end+1) = skelPix / hullPix;

            if ~isempty(juncXY)
                inTile = juncXY(:,1) >= cc(1) & juncXY(:,1) <= cc(end) & ...
                         juncXY(:,2) >= rr(1) & juncXY(:,2) <= rr(end);
                jdens(end+1) = sum(inTile) / hullPix; %#ok<AGROW>
            end
        end
    end

    nUsed = numel(dens);
    if nUsed >= 2
        cvDens = std(dens) / max(mean(dens), eps);
        if numel(jdens) >= 2
            cvJunc = std(jdens) / max(mean(jdens), eps);
        end
    end
end

% ----------------------------------------------------------------------
function [cvDens, profile] = radialHeterogeneity(skel_cap, im_hull, apexXY, nBins)
% Vessel density in equal-width distance bands from the optic-nerve apex.
% Retinal capillary density varies systematically with eccentricity, so this
% is the anatomically meaningful axis - but it needs the apex location.

    cvDens = NaN; profile = nan(1, nBins);
    if isempty(apexXY), return; end

    [rr, cc] = size(im_hull);
    [X, Y]   = meshgrid(1:cc, 1:rr);
    Dap      = hypot(X - apexXY(1), Y - apexXY(2));

    dHull = Dap(im_hull);
    if isempty(dHull), return; end
    edges = linspace(min(dHull), max(dHull)+eps, nBins+1);

    for b = 1:nBins
        band = im_hull & Dap >= edges(b) & Dap < edges(b+1);
        hp   = nnz(band);
        if hp > 0
           profile(b) = nnz(skel_cap & band) / hp;
        end
    end

    good = profile(isfinite(profile));
    if numel(good) >= 2
        cvDens = std(good) / max(mean(good), eps);
    end
end

% ----------------------------------------------------------------------
function s = emptyStats()
% Single source of truth for the scalar field list, so the batch script can
% build its results table without hand-listing every name.
    f = {'nJunctions','nEndpoints','nBranches','nComponents', ...
         'leafletArea_mm2','capSkelLen_mm', ...
         'junctionsPerArea_mm2','junctionsPerLength_mm','endpointJunctionRatio', ...
         'holeCount','holesPerArea_mm2','holeCountMask', ...
         'meanBranchLen_um','sdBranchLen_um','medianBranchLen_um', ...
         'iqrBranchLen_um','cvBranchLen', ...
         'meanTortuosity','sdTortuosity','medianTortuosity', ...
         'nHoles','meanHoleArea_um2','medianHoleArea_um2','sdHoleArea_um2', ...
         'maxHoleArea_um2','cvHoleArea','avascularFrac', ...
         'meanIntercapDist_um','p95IntercapDist_um','maxIntercapDist_um', ...
         'meanCapRadius_um','sdCapRadius_um', ...
         'gridCV_density','gridCV_junctionDensity','nGridTiles','radialCV_density', ...
         'largeSkelPixels','capSkelPixels','radialDens_b1','radialDens_b2','radialDens_b3', ...
         'radialDens_b4','radialDens_b5','radialDens_b6'};
    for k = 1:numel(f), s.(f{k}) = NaN; end
end