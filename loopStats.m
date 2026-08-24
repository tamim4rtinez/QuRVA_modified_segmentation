function [Tret, Tbin] = loopStats(masterFolder, params)
% loopStats  Capillary mesh loop metrics, gap-closed and size-filtered.
%
%   [Tret, Tbin] = loopStats(masterFolder)
%   loopStats(masterFolder, struct('gapPix', 2, 'minLoopUM2', 100))
%
% Replaces the original nHoles / holeCount / holeCountMask / HoleArea
% columns, which were computed on an unclosed skeleton and an unfiltered
% fill - roughly 80% of those "loops" were single-pixel diagonal
% staircase artifacts, and genuine loops left open by one-pixel gaps
% were missed entirely.
%
% HERE: the capillary skeleton is closed by gapPix before filling, so
% near-miss loops close; enclosed regions below minLoopUM2 are dropped
% as artifacts. Loop areas are measured to the skeleton CENTRELINES,
% so they run larger than mask-based hole areas by about one vessel
% width - not comparable to the old columns.
%
% Writes loop_stats_retina.csv (one row per retina) and
%        loop_stats_bins.csv   (one row per retina per radial bin)

    if nargin < 2, params = struct(); end
    p.basePxUM       = getf(params,'basePxUM',       3/1.4491);
    p.radiusThreshUM = getf(params,'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params,'minBranchPix',   8);
    p.minLoopUM2     = getf(params,'minLoopUM2',     100);
    p.gapPix         = getf(params,'gapPix',         2);
    p.nBins          = getf(params,'nBins',          6);

    files = dir(fullfile(masterFolder,'VasculatureNumbers','*.mat'));
    assert(~isempty(files), 'No .mat files found');

    Rrows = []; Brows = [];

    for ii = 1:numel(files)
        f = files(ii).name; base = string(erase(f,'.mat'));
        try
            S = load(fullfile(masterFolder,'VasculatureNumbers',f), ...
                     'smoothVessels','vesselSkelMask');
            M = load(fullfile(masterFolder,'Masks',f), 'thisMask');
            C = load(fullfile(masterFolder,'ONCenter',f),'thisONCenter');

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

            D  = bwdist(~im_clean);
            sr = zeros(size(D)); sr(im_skel) = D(im_skel);
            skel_cap = bwareaopen(im_skel & sr > 0 & ...
                                  sr <= p.radiusThreshUM/pxUM, p.minBranchPix, 8);

            % --- close one-pixel gaps, then re-thin to a skeleton ---
            skel_fix = bwmorph(imclose(skel_cap, strel('disk', p.gapPix)), ...
                               'thin', Inf);

            allLoops = imfill(skel_fix,'holes') & ~skel_fix;
            minPix   = max(round(p.minLoopUM2 / pxUM^2), 2);

            cc   = bwconncomp(allLoops, 8);
            Apix = cellfun(@numel, cc.PixelIdxList)';
            keep = Apix >= minPix;
            nRej = nnz(~keep);

            props = regionprops(cc, 'Centroid');
            allC  = cat(1, props.Centroid);
            loopA = Apix(keep) * pxUM^2;
            loopC = allC(keep, :);

            hullArea_mm2 = nnz(im_hull) * pxUM^2 / 1e6;

            % ---------- whole retina ----------
            r.fileName            = base;
            r.strain              = string(extractBefore(char(base),'_'));
            r.pxUM                = pxUM;
            r.hullArea_mm2        = hullArea_mm2;
            r.nLoops              = numel(loopA);
            r.nRejected           = nRej;
            r.fracKept            = numel(loopA) / max(numel(Apix),1);
            r.loopDensity_mm2     = numel(loopA) / max(hullArea_mm2, eps);
            r.meanLoopArea_um2    = safe(@mean,   loopA);
            r.medianLoopArea_um2  = safe(@median, loopA);
            r.sdLoopArea_um2      = safe(@std,    loopA);
            r.iqrLoopArea_um2     = safe(@iqr,    loopA);
            r.p95LoopArea_um2     = safe(@(x) prctile(x,95), loopA);
            r.maxLoopArea_um2     = safe(@max,    loopA);
            r.cvLoopArea          = safe(@std, loopA) / max(safe(@mean, loopA), eps);
            r.totalLoopArea_mm2   = sum(loopA) / 1e6;
            r.loopAreaFrac        = (sum(loopA)/1e6) / max(hullArea_mm2, eps);
            Rrows = [Rrows; r]; %#ok<AGROW>

            % ---------- per radial bin ----------
            [nr, nc] = size(im_hull);
            [X,Y]    = meshgrid(1:nc,1:nr);
            Dap      = hypot(X-apexXY(1), Y-apexXY(2));
            dH       = Dap(im_hull);
            edges    = linspace(min(dH), max(dH)+eps, p.nBins+1);
            loopD    = hypot(loopC(:,1)-apexXY(1), loopC(:,2)-apexXY(2));

            for b = 1:p.nBins
                band  = im_hull & Dap >= edges(b) & Dap < edges(b+1);
                a_mm2 = nnz(band) * pxUM^2 / 1e6;
                inB   = loopD >= edges(b) & loopD < edges(b+1);
                la    = loopA(inB);

                q.fileName           = base;
                q.strain             = r.strain;
                q.bin                = b;
                q.inner_um           = edges(b)*pxUM;
                q.outer_um           = edges(b+1)*pxUM;
                q.tissueArea_mm2     = a_mm2;
                q.nLoops             = numel(la);
                q.loopDensity_mm2    = numel(la) / max(a_mm2, eps);
                q.meanLoopArea_um2   = safe(@mean,   la);
                q.medianLoopArea_um2 = safe(@median, la);
                q.iqrLoopArea_um2    = safe(@iqr,    la);
                q.p95LoopArea_um2    = safe(@(x) prctile(x,95), la);
                q.cvLoopArea         = safe(@std, la) / max(safe(@mean, la), eps);
                q.loopAreaFrac       = (sum(la)/1e6) / max(a_mm2, eps);
                Brows = [Brows; q]; %#ok<AGROW>
            end

            fprintf('[%2d/%2d] %-40s px %.3f  %5d loops (%4.0f%% kept)  med %5.0f um2\n', ...
                ii, numel(files), base, pxUM, r.nLoops, 100*r.fracKept, ...
                r.medianLoopArea_um2);

        catch ME
            fprintf('[%2d/%2d] %-40s FAILED: %s\n', ii, numel(files), base, ME.message);
        end
    end

    Tret = struct2table(Rrows);
    Tbin = struct2table(Brows);
    writetable(Tret, fullfile(masterFolder,'loop_stats_retina.csv'));
    writetable(Tbin, fullfile(masterFolder,'loop_stats_bins.csv'));

    fprintf('\nWrote loop_stats_retina.csv (%d rows) and loop_stats_bins.csv (%d rows)\n', ...
            height(Tret), height(Tbin));
    fprintf('\n--- fraction of enclosed regions kept, by strain ---\n');
    disp(groupsummary(Tret, 'strain', {'mean','min'}, 'fracKept'));
end

function v = safe(fn,x), if isempty(x), v = NaN; else, v = fn(double(x)); end, end
function v = getf(s,f,d), if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end, end