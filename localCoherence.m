function [T, maps] = localCoherence(masterFolder, fileName, params)
% localCoherence  Local vessel alignment in small windows, summarised per bin.
%
%   [T, maps] = localCoherence(masterFolder, 'ACI_1_....tif.mat')
%   localCoherence(..., struct('showQC', true))
%
% WHY THIS AND NOT radialAlignment:
%   radialAlignment averages over a whole annulus. If a retina has patches
%   of parallel vessels pointing in DIFFERENT directions, those patches
%   cancel and the annulus reads zero - the same as a uniform mesh.
%   Here alignment is computed in small windows first. A patchy-parallel
%   retina then shows a TAIL of high-coherence windows even though the
%   mean over the bin is unremarkable.
%
% METHOD
%   Tile the image with overlapping windows of side tileUM.
%   In each window, over the branches whose centre falls inside it:
%       localCoherence = |sum(w * exp(2i*theta))| / sum(w)
%   theta = branch principal-axis angle, w = branch length in pixels.
%   Angle doubling treats a branch as an undirected line.
%   0 = window has no common direction, 1 = every branch in the window
%   is parallel. Windows with fewer than minBranches are skipped.
%
% PER-BIN SUMMARIES (a window belongs to the bin its centre falls in)
%   meanCoh     average window alignment
%   p90Coh      90th percentile - sensitive to a minority of aligned patches
%   fracAligned fraction of windows above cohThresh  <- the patchiness measure
%   nTiles      windows contributing (under ~30 and the bin is unreliable

    if nargin < 3, params = struct(); end
    p.basePxUM       = getf(params,'basePxUM',       3/1.4491);
    p.nBins          = getf(params,'nBins',          6);
    p.radiusThreshUM = getf(params,'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params,'minBranchPix',   8);
    p.minLenUM       = getf(params,'minLenUM',       20);
    p.tileUM         = getf(params,'tileUM',         250);
    p.minBranches    = getf(params,'minBranches',    6);
    p.cohThresh      = getf(params,'cohThresh',      0.60);
    p.showQC         = getf(params,'showQC',         true);

    S = load(fullfile(masterFolder,'VasculatureNumbers',fileName), ...
             'smoothVessels','vesselSkelMask');
    M = load(fullfile(masterFolder,'Masks',fileName),   'thisMask');
    C = load(fullfile(masterFolder,'ONCenter',fileName),'thisONCenter');

    im_clean = logical(S.smoothVessels);
    im_skel  = logical(S.vesselSkelMask);
    im_hull  = logical(M.thisMask);


    sf = 1;
    if ~isequal(size(im_hull), size(im_clean))
        sf = size(im_hull,1)/size(im_clean,1);
        im_hull = imresize(im_hull, size(im_clean), 'nearest');
    end
    apexXY = C.thisONCenter(:).' / sf;
    pxUM   = p.basePxUM * sf;

    D = bwdist(~im_clean);
    sr = zeros(size(D)); sr(im_skel) = D(im_skel);
    skel_cap = bwareaopen(im_skel & sr > 0 & sr <= p.radiusThreshUM/pxUM, ...
                          p.minBranchPix, 8);

    nbr    = conv2(double(skel_cap), [1 1 1;1 0 1;1 1 1], 'same');
    L      = bwlabel(skel_cap & ~(skel_cap & nbr >= 3), 8);
    st     = regionprops(L, 'PixelList','Centroid');

    minPix = max(round(p.minLenUM/pxUM), 5);

    cen = []; th = []; wt = [];
    for k = 1:numel(st)
        P = st(k).PixelList;
        if size(P,1) < minPix, continue; end
        Cv = cov(double(P));
        [V,Dg] = eig(Cv); [~,i2] = max(diag(Dg));
        v = V(:,i2)/norm(V(:,i2));
        cen(end+1,:) = st(k).Centroid;      %#ok<AGROW>
        th(end+1,1)  = atan2(v(2), v(1));   %#ok<AGROW>
        wt(end+1,1)  = size(P,1);           %#ok<AGROW>
    end

    % ---- tiles, 50% overlap ----
    [nr, nc] = size(im_hull);
    tile = max(round(p.tileUM/pxUM), 8);
    step = round(tile/2);

    tx = []; ty = []; coh = []; nb = [];
    for y0 = 1:step:(nr-tile+1)
        for x0 = 1:step:(nc-tile+1)
            in = cen(:,1) >= x0 & cen(:,1) < x0+tile & ...
                 cen(:,2) >= y0 & cen(:,2) < y0+tile;
            if nnz(in) < p.minBranches, continue; end
            w = wt(in);
            c = abs(sum(w .* exp(2i*th(in)))) / sum(w);
            nb(end+1,1) = nnz(in);   %#ok<AGROW>
            tx(end+1,1)  = x0 + tile/2;  %#ok<AGROW>
            ty(end+1,1)  = y0 + tile/2;  %#ok<AGROW>
            coh(end+1,1) = c;            %#ok<AGROW>
        end
    end

    % ---- assign tiles to radial bins ----
    [X,Y]  = meshgrid(1:nc,1:nr);
    Dap    = hypot(X-apexXY(1), Y-apexXY(2));
    dHull  = Dap(im_hull);
    edges  = linspace(min(dHull), max(dHull)+eps, p.nBins+1);
    tileD  = hypot(tx-apexXY(1), ty-apexXY(2));

    base = string(erase(fileName,'.mat'));
    rows = [];
    fprintf('\n%s   px %.3f um   tile %d px (%.0f um)   %d tiles total\n', ...
            base, pxUM, tile, tile*pxUM, numel(coh));
    for b = 1:p.nBins
        in = tileD >= edges(b) & tileD < edges(b+1);
        c  = coh(in);
        r.fileName    = base;
        r.bin         = b;
        r.inner_um    = edges(b)*pxUM;
        r.outer_um    = edges(b+1)*pxUM;
        r.nTiles      = numel(c);
        r.meanCoh     = safe(@mean, c);
        r.p90Coh      = safe(@(x) prctile(x,90), c);
        r.fracAligned = safe(@(x) mean(x > p.cohThresh), c);
        rows = [rows; r]; %#ok<AGROW>
        fprintf('  bin %d  %6.0f-%6.0f um  n %4d  mean %.3f  p90 %.3f  frac>%.2f %.3f\n', ...
            b, r.inner_um, r.outer_um, r.nTiles, r.meanCoh, r.p90Coh, ...
            p.cohThresh, r.fracAligned);
    end
    T = struct2table(rows);
        maps = struct('tx',tx,'ty',ty,'coh',coh,'nbr',nb, ...
                  'edges',edges,'pxUM',pxUM,'apex',apexXY);

    if p.showQC
        figure('Color','w','Name','local coherence map');
        imshow(0.15*double(cat(3,skel_cap,skel_cap,skel_cap))+0.04); hold on
        scatter(tx, ty, (tile/6)^2, coh, 'filled', 'MarkerFaceAlpha', 0.55);
        colormap(gca, hot); caxis([0 1]); cb = colorbar; cb.Color = 'k';
        cb.Label.String = 'local coherence (1 = parallel)';
        a = linspace(0,2*pi,400);
        for b = 2:p.nBins
            plot(apexXY(1)+edges(b)*cos(a), apexXY(2)+edges(b)*sin(a), 'w-','LineWidth',0.8);
        end
        plot(apexXY(1), apexXY(2), 'w*','MarkerSize',12,'LineWidth',1.2);
        title(sprintf('%s  bright = locally parallel', strrep(base,'_','\_')));

        figure('Color','w','Name','coherence distribution');
        histogram(coh, 0:0.05:1, 'Normalization','probability', ...
                  'FaceColor',[0.35 0.45 0.7]);
        xlabel('local coherence'); ylabel('fraction of windows');
        title(sprintf('%s  (all bins pooled)', strrep(base,'_','\_')));
    end
end

function v = safe(fn,x), if isempty(x), v = NaN; else, v = fn(double(x)); end, end
function v = getf(s,f,d), if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end, end