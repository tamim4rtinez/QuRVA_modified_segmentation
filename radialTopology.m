function T = radialTopology(masterFolder, params)
% radialTopology  Capillary network topology in radial bins from the ONH.
%
%   T = radialTopology(masterFolder)
%   T = radialTopology(masterFolder, params)
%
% Everything is computed on the CAPILLARY SKELETON (skel_cap), never the
% vessel mask - except intercapillary distance, which by definition needs
% the mask (it measures tissue distance to the nearest vessel).
%
% Bins are equal-width fractions of each retina's own min-to-max distance
% from the optic nerve head, matching capillaryMetrics/radialHeterogeneity.
%
% params (all optional):
%   basePxUM         um per pixel of the segmented TIFF.  Default 3/1.4491
%   nBins            radial bins.                          Default 6
%   radiusThreshUM   capillary/large-vessel cutoff, RADIUS. Default 17.1
%                    (= 4.5 px at the rat effective pixel size, i.e. the
%                     same split the existing rat table already used)
%   minBranchPix     spur removal.                          Default 8
%
% Output: one row per retina per bin
%   icdMean_um, icdP95_um        intercapillary distance (mask-based)
%   loopDensity_mm2              Euler loops per mm2 of tissue
%   medianLoopArea_um2           median area of enclosed skeleton loops
%   nLoops                       number of enclosed loops measured
%   meanBranchLen_um, medianBranchLen_um, iqrBranchLen_um, cvBranchLen
%   nBranches
%   endpointJunctionRatio, nEndpoints, nJunctions
%   capLen_mm, tissueArea_mm2, coverage

    if nargin < 2, params = struct(); end
    p.basePxUM       = getf(params, 'basePxUM',       3/1.4491);
    p.nBins          = getf(params, 'nBins',          6);
    p.radiusThreshUM = getf(params, 'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params, 'minBranchPix',   8);

    vascDir = fullfile(masterFolder, 'VasculatureNumbers');
    maskDir = fullfile(masterFolder, 'Masks');
    ocDir   = fullfile(masterFolder, 'ONCenter');
    assert(exist(vascDir,'dir')==7, 'No VasculatureNumbers in %s', masterFolder);
    assert(exist(ocDir,'dir')==7, ...
        ['No ONCenter folder - radial binning needs the optic nerve head ' ...
         'position from QuRVA.']);

    files = dir(fullfile(vascDir, '*.mat'));
    assert(~isempty(files), 'No .mat files in %s', vascDir);
    fprintf('%d retinas, %d radial bins\n', numel(files), p.nBins);

    rows = [];

    for ii = 1:numel(files)
        base = erase(files(ii).name, '.mat');

        try
            S = load(fullfile(vascDir, files(ii).name), ...
                     'smoothVessels','vesselSkelMask');
            M = load(fullfile(maskDir, files(ii).name), 'thisMask');
            C = load(fullfile(ocDir,   files(ii).name), 'thisONCenter');

            im_clean = logical(S.smoothVessels);
            im_skel  = logical(S.vesselSkelMask);
            im_hull  = logical(M.thisMask);

            sf = 1;
            if ~isequal(size(im_hull), size(im_clean))
                sf = size(im_hull,1) / size(im_clean,1);
                im_hull = imresize(im_hull, size(im_clean), 'nearest');
            end
            apexXY = C.thisONCenter(:).' / sf;
            pxUM   = p.basePxUM * sf;

            % ---------- capillary skeleton ----------
            D = bwdist(~im_clean);
            skelRadius = zeros(size(D));
            skelRadius(im_skel) = D(im_skel);
            threshPix = p.radiusThreshUM / pxUM;
            skel_cap  = im_skel & skelRadius > 0 & skelRadius <= threshPix;
            skel_cap  = bwareaopen(skel_cap, p.minBranchPix, 8);

            nbr    = conv2(double(skel_cap), [1 1 1;1 0 1;1 1 1], 'same');
            juncPx = skel_cap & nbr >= 3;
            endPx  = skel_cap & nbr == 1;

            % branch segments, labelled once for the whole retina
            branches = skel_cap & ~juncPx;
            L = bwlabel(branches, 8);
            [brLen_px, brCentroid] = branchLengths(L);
            brLen_um = brLen_px * pxUM;

            % enclosed loops of the skeleton, as objects
            % for one rat image, with skel_cap already built

            loopMask = imfill(skel_cap, 'holes') & ~skel_cap;   % 4-conn, as originally
            minLoopPix = round(100 / pxUM^2);
            loopMask = bwareaopen(loopMask, max(minLoopPix, 2), 8);
            ccLoop   = bwconncomp(loopMask, 8);
            if ccLoop.NumObjects > 0
                lp = regionprops(ccLoop, 'Area','Centroid');
                loopArea_um2 = [lp.Area]' * pxUM^2;
                loopXY       = cat(1, lp.Centroid);      % [x y]
            else
                loopArea_um2 = zeros(0,1); loopXY = zeros(0,2);
            end

            % ---------- radial bands ----------
            [nr, nc] = size(im_hull);
            [X, Y]   = meshgrid(1:nc, 1:nr);
            Dap      = hypot(X - apexXY(1), Y - apexXY(2));
            dHull    = Dap(im_hull);
            edges    = linspace(min(dHull), max(dHull)+eps, p.nBins+1);

            Dout     = bwdist(im_clean);       % tissue -> nearest vessel
            inTissue = im_hull & ~im_clean;

              

            % distance of each branch / loop centroid from the apex
            brD   = hypot(brCentroid(:,1) - apexXY(1), brCentroid(:,2) - apexXY(2));
            loopD = hypot(loopXY(:,1)     - apexXY(1), loopXY(:,2)     - apexXY(2));

            for b = 1:p.nBins
                band  = im_hull & Dap >= edges(b) & Dap < edges(b+1);
                ring  = Dap >= edges(b) & Dap < edges(b+1);

                tissuePix = nnz(band);
                area_mm2  = tissuePix * pxUM^2 / 1e6;

                icd = Dout(band & ~im_clean) * pxUM;

                nJ = nnz(juncPx & band);
                nE = nnz(endPx  & band);

                inBr = brD >= edges(b) & brD < edges(b+1);
                bl   = brLen_um(inBr);

                inLp = loopD >= edges(b) & loopD < edges(b+1);
                la   = loopArea_um2(inLp);

                nLoopsBand = numel(la);

                r.fileName  = string(base);
                r.bin       = b;
                r.pxUM      = pxUM;
                r.tissueArea_mm2 = area_mm2;
                r.coverage  = tissuePix / max(nnz(ring), 1);

                r.icdMean_um = safe(@mean, icd);
                r.icdP95_um  = safe(@(x) prctile(x,95), icd);

                r.nLoops            = nLoopsBand;
                r.loopDensity_mm2   = nLoopsBand / max(area_mm2, eps);
                r.medianLoopArea_um2 = safe(@median, la);

                r.nBranches          = numel(bl);
                r.meanBranchLen_um   = safe(@mean,   bl);
                r.medianBranchLen_um = safe(@median, bl);
                r.iqrBranchLen_um    = safe(@iqr,    bl);
                r.cvBranchLen        = safe(@std, bl) / max(safe(@mean, bl), eps);

                r.nJunctions = nJ;
                r.nEndpoints = nE;
                r.endpointJunctionRatio = nE / max(nJ, 1);

                r.capLen_mm = sum(bl) / 1e3;
                ro = radialOrientation(skel_cap, apexXY, band, pxUM, 20);
                r.radialAlignment     = ro.radialAlignment;
                r.globalCoherence     = ro.globalCoherence;
                r.meanRadialAngle_deg = ro.meanRadialAngle_deg;
                r.nOriented           = ro.nUsed;

                rows = [rows; r]; %#ok<AGROW>
            end

            fprintf('[%3d/%3d] %-42s px %.3f um  loops/bin:', ...
                ii, numel(files), base, pxUM);
                fprintf(' %4d', [rows(end-p.nBins+1:end).nLoops]);
            fprintf('\n');

        catch ME
            fprintf('[%3d/%3d] %-42s FAILED: %s\n', ...
                ii, numel(files), base, ME.message);
        end
    end

    T = struct2table(rows);
    out = fullfile(masterFolder, 'radial_topology.csv');
    writetable(T, out);
    fprintf('\nWrote %s\n', out);
end

% ======================================================================
function [len, cent] = branchLengths(L)
% Path length of each labelled branch, diagonal steps as sqrt(2),
% plus each branch's centroid so it can be assigned to a radial band.
    n = max(L(:));
    len  = zeros(n,1);
    offs = {[0 1],[1 0],[1 1],[1 -1]};
    step = [1, 1, sqrt(2), sqrt(2)];
    for k = 1:4
        Ls = circshift(L, offs{k});
        if offs{k}(1) ~= 0,  Ls(1,:)   = 0; end
        if offs{k}(2) ==  1, Ls(:,1)   = 0; end
        if offs{k}(2) == -1, Ls(:,end) = 0; end
        same = (L > 0) & (L == Ls);
        if any(same(:))
            len = len + accumarray(L(same), step(k), [n 1]);
        end
    end
    pix = accumarray(L(L>0), 1, [n 1]);
    len(len == 0 & pix > 0) = 1;

    if n > 0
        st = regionprops(L, 'Centroid');
        cent = cat(1, st.Centroid);        % [x y]
    else
        cent = zeros(0,2);
    end
end

function v = safe(fn, x)
    if isempty(x), v = NaN; else, v = fn(double(x)); end
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end