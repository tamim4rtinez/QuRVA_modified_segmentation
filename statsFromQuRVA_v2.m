function T = statsFromQuRVA_v2(masterFolder, capP)
% statsFromQuRVA  Whole-retina capillary metrics from QuRVA's segmentations.
%
%   T = statsFromQuRVA(masterFolder, capP)
%
% Reads the .mat files QuRVA writes to VasculatureNumbers/ - these hold the
% RAW BINARY masks, so nothing is re-segmented and nothing is recovered from
% the display overlays in VasculatureImages/:
%     smoothVessels    vessel mask
%     vesselSkelMask   skeleton
%     brchPts/endPts   QuRVA's own branch and end points, kept for comparison
% Retinal area comes from Masks/thisMask
% 
% 
% , which is tighter than a convex hull.
%
% Point masterFolder at the IMAGE folder you gave processFolder, not at a
% subfolder.
%
% ---------------------------------------------------------------------
% REQUIRED
%   capP.pixSzeUM        microns per pixel of the TIFF you gave
%                        processFolder. Every length and area output
%                        depends on it. There is no safe default.
%
% OPTIONAL
%   capP.radiusThreshUM  capillary / large-vessel cutoff, RADIUS in um.
%                        Preferred over radiusThreshPix: it means the same
%                        physical vessel size regardless of image scale.
%   capP.radiusThreshPix used only if radiusThreshUM is not given. Note this
%                        is a different physical size in every image whose
%                        resolution differs.       default 4.5
%   capP.minBranchPix    spur removal.             default 8
%   capP.gridN           grid tiles per side.      default 8
%   capP.gapPix          skeleton gap closing before loop detection.
%                        default 2. See LOOPS below - worth tuning per
%                        dataset.
%   capP.minLoopUM2      smallest enclosed region counted as a real mesh
%                        loop.                     default 100
%   capP.qcScale         QC image scale, 1 = full size.  default 0.5
%   capP.debug           true rethrows errors instead of logging them.
%
% ---------------------------------------------------------------------
% LOOPS
% Mesh loops are computed here, not in capillaryMetrics, because the naive
% version is badly wrong in two directions at once:
%
%   - A skeleton drawn on a pixel grid turns diagonals into staircases, and
%     every notch registers as a tiny enclosed region. Unfiltered, roughly
%     80% of detected "loops" are these artifacts.
%   - Skeletonisation sometimes leaves a one-pixel break in a vessel, which
%     lets the fill leak into what should be a closed loop, so genuine mesh
%     spaces go uncounted.
%
% So the capillary skeleton is closed by gapPix and re-thinned before
% filling, and enclosed regions below minLoopUM2 are dropped.
%
% TUNING gapPix: increase it until the loop count plateaus. If the count
% keeps climbing, adjacent loops are being merged and you have gone too far.
% TUNING minLoopUM2: plot a histogram of log10 loop area - the artifact
% population and the real population separate cleanly, and the valley
% between them is the threshold.
%
% Loop areas are measured to skeleton CENTRELINES, so they run larger than
% mask-based hole areas by roughly one vessel width.
%
% The old holeCount / nHoles / holesPerArea_mm2 / *HoleArea_um2 columns from
% capillaryMetrics are removed from the output: they were the unclosed,
% unfiltered version, and several were computed on the full vessel mask
% rather than the capillary skeleton.
% ---------------------------------------------------------------------

    if nargin < 2, capP = struct(); end

    if ~isfield(capP,'pixSzeUM') || isempty(capP.pixSzeUM)
        error(['capP.pixSzeUM is required.\n' ...
               'It is the microns per pixel of the TIFF you gave ' ...
               'processFolder - not of the segmented image, which QuRVA ' ...
               'rescales internally.\nEvery length, area and density ' ...
               'output is wrong without it, silently.']);
    end

    p = capP;
    if ~isfield(p,'radiusThreshPix'), p.radiusThreshPix = 4.5;  end
    if ~isfield(p,'minBranchPix'),    p.minBranchPix    = 8;    end
    if ~isfield(p,'gridN'),           p.gridN           = 8;    end
    if ~isfield(p,'gapPix'),          p.gapPix          = 2;    end
    if ~isfield(p,'minLoopUM2'),      p.minLoopUM2      = 100;  end
    if ~isfield(p,'qcScale'),         p.qcScale         = 0.5;  end
    if ~isfield(p,'debug'),           p.debug           = false;end

    vascDir = fullfile(masterFolder, 'VasculatureNumbers');
    maskDir = fullfile(masterFolder, 'Masks');
    qcDir   = fullfile(masterFolder, 'CapillaryQC');

    if ~exist(vascDir, 'dir')
        error(['No VasculatureNumbers folder in %s\n' ...
               'Point this at the image folder you gave processFolder.'], ...
               masterFolder);
    end
    if ~exist(qcDir,'dir'), mkdir(qcDir); end

    files = dir(fullfile(vascDir, '*.mat'));
    if isempty(files), error('No .mat files in %s', vascDir); end
    fprintf('Found %d processed images\n', numel(files));

    % Fields capillaryMetrics returns, minus the hole family we replace
    capTemplate = capillaryMetrics();
    capFields   = fieldnames(capTemplate);
       dropFields  = {'holeCount','nHoles','holesPerArea_mm2','holeCountMask', ...
                   'meanHoleArea_um2','medianHoleArea_um2','sdHoleArea_um2', ...
                   'maxHoleArea_um2','cvHoleArea', ...
                   'gridCV_density','gridCV_junctionDensity','nGridTiles', ...
                   'radialCV_density', 'radialDens_b1',	'radialDens_b2', 'radialDens_b3', 'radialDens_b4', 'radialDens_b5', 'radialDens_b6'};
    capFields   = setdiff(capFields, dropFields, 'stable');

    loopFields = {'nLoops','nRejectedLoops','fracLoopsKept', ...
                  'loopDensity_mm2','meanLoopArea_um2','medianLoopArea_um2', ...
                  'sdLoopArea_um2','iqrLoopArea_um2','p95LoopArea_um2', ...
                  'maxLoopArea_um2','cvLoopArea','totalLoopArea_mm2', ...
                  'loopAreaFrac'};

        core = struct('fileName','', 'ok',false, 'msg','', ...
                  'pixSzeUM_used',NaN, 'radiusThreshUM_used',NaN, ...
                  'maskPixels',NaN, 'skelPixels',NaN, 'allSkelLen_mm',NaN, ...
                  'maskFrac',NaN, ...
                  'qurvaBranchPts',NaN, 'qurvaEndPts',NaN, ...
                  'qurvaAvascArea',NaN, 'scaleFactor',NaN);
    for f = 1:numel(capFields),  core.(capFields{f})  = NaN; end
    for f = 1:numel(loopFields), core.(loopFields{f}) = NaN; end
    results = repmat(core, numel(files), 1);

    for ii = 1:numel(files)

        base = erase(files(ii).name, '.mat');   % keeps the original extension
        results(ii).fileName = base;

        try
            S = load(fullfile(vascDir, files(ii).name), ...
                     'smoothVessels','vesselSkelMask','brchPts','endPts','aVascZone');

            im_clean = logical(S.smoothVessels);
            im_skel  = logical(S.vesselSkelMask);

            % ---------- retinal mask and scale ----------
            % QuRVA rescales the image (resetScale) before segmenting, so the
            % saved retinal mask may be at a different resolution than the
            % vessel mask. Match them, and record the factor - it also scales
            % the true micron-per-pixel value.
            hullMask    = [];
            scaleFactor = 1;
            mp = fullfile(maskDir, files(ii).name);
            if exist(mp, 'file')
                M = load(mp, 'thisMask');
                hullMask = logical(M.thisMask);
                if ~isequal(size(hullMask), size(im_clean))
                    scaleFactor = size(hullMask,1) / size(im_clean,1);
                    hullMask = imresize(hullMask, size(im_clean), 'nearest');
                end
            end
            results(ii).scaleFactor = scaleFactor;

            pxUM = p.pixSzeUM * scaleFactor;     % um per pixel of im_clean
            results(ii).pixSzeUM_used = pxUM;

            % radius threshold: physical if given, else pixels
            if isfield(p,'radiusThreshUM') && ~isempty(p.radiusThreshUM)
                threshPix = p.radiusThreshUM / pxUM;
            else
                threshPix = p.radiusThreshPix;
            end
            results(ii).radiusThreshUM_used = threshPix * pxUM;

            % ---------- metrics ----------
            capPi = p;
            capPi.hull            = hullMask;    % [] falls back to bwconvhull
            capPi.pixSzeUM        = pxUM;
            capPi.radiusThreshPix = threshPix;

            [cap, capQC] = capillaryMetrics(im_clean, im_skel, capPi);

            for f = 1:numel(capFields)
                results(ii).(capFields{f}) = cap.(capFields{f});
            end

            % ---------- mesh loops (see LOOPS in the header) ----------
            lp = loopMetrics(im_clean, im_skel, hullMask, pxUM, threshPix, p);
            for f = 1:numel(loopFields)
                results(ii).(loopFields{f}) = lp.(loopFields{f});
            end

            % ---------- bookkeeping ----------
            results(ii).ok             = true;
            results(ii).maskPixels     = nnz(im_clean);
            results(ii).skelPixels     = nnz(im_skel);
                        results(ii).allSkelLen_mm  = skelLength(im_skel) * pxUM / 1e3;
            if ~isempty(hullMask)
                results(ii).maskFrac = nnz(im_clean) / max(nnz(hullMask),1);
            end
            results(ii).qurvaBranchPts = nnz(S.brchPts);
            results(ii).qurvaEndPts    = nnz(S.endPts);
            if isfield(S,'aVascZone'), results(ii).qurvaAvascArea = nnz(S.aVascZone); end

            if p.qcScale < 1, capQC = imresize(capQC, p.qcScale); end
            imwrite(capQC, fullfile(qcDir, [base '_cap.png']));

            fprintf(['[%3d/%3d] %-45s ok  px %.3f um  junc %5d  ' ...
                     'loops %5d (%3.0f%% kept)  meanBr %6.1f\n'], ...
                ii, numel(files), base, pxUM, cap.nJunctions, ...
                lp.nLoops, 100*lp.fracLoopsKept, cap.meanBranchLen_um);

        catch ME
            if p.debug, rethrow(ME); end
            results(ii).ok  = false;
            results(ii).msg = ME.message;
            fprintf('[%3d/%3d] %-45s FAILED: %s\n', ii, numel(files), base, ME.message);
        end
    end

    T = struct2table(results);
    csvPath = fullfile(masterFolder, 'qurva_capillary_stats.csv');
    writetable(T, csvPath);

    ok = [results.ok];
    fprintf('\n%d/%d succeeded.\nStats: %s\nQC: %s\n', ...
        sum(ok), numel(files), csvPath, qcDir);

    if any(ok)
        sf = [results(ok).scaleFactor];
        if numel(unique(round(sf,3))) > 1
            warning(['scaleFactor is not constant across images (%s).\n' ...
                     'That means the source images were not all at the same ' ...
                     'resolution - most often because they were exported from ' ...
                     'different .czi pyramid series. A single pixSzeUM cannot ' ...
                     'be correct for all of them.'], ...
                     mat2str(unique(round(sf,3))));
        end
        fk = [results(ok).fracLoopsKept];
        if median(fk) < 0.5
            warning(['Over half of enclosed regions are being rejected as too ' ...
                     'small (median kept %.0f%%). Check the loop-area histogram ' ...
                     'and consider raising gapPix or lowering minLoopUM2.'], ...
                     100*median(fk));
        end
    end
end

% ======================================================================
function lp = loopMetrics(im_clean, im_skel, hullMask, pxUM, threshPix, p)
% Mesh loops of the CAPILLARY skeleton, gap-closed and size-filtered.

    D  = bwdist(~im_clean);
    sr = zeros(size(D)); sr(im_skel) = D(im_skel);
    skel_cap = bwareaopen(im_skel & sr > 0 & sr <= threshPix, p.minBranchPix, 8);

    % Close one-pixel breaks so near-miss loops close, then re-thin so the
    % result is still a skeleton rather than a thickened mask.
    skel_fix = bwmorph(imclose(skel_cap, strel('disk', p.gapPix)), 'thin', Inf);

    % 4-connected background fill: enclosed regions of the skeleton
    allLoops = imfill(skel_fix, 'holes') & ~skel_fix;

    minPix = max(round(p.minLoopUM2 / pxUM^2), 2);
    cc     = bwconncomp(allLoops, 8);
    Apix   = cellfun(@numel, cc.PixelIdxList)';
    keep   = Apix >= minPix;

    loopA = Apix(keep) * pxUM^2;                       % um^2

    if isempty(hullMask), hullMask = bwconvhull(im_clean); end
    hullArea_mm2 = nnz(hullMask) * pxUM^2 / 1e6;

    lp.nLoops             = numel(loopA);
    lp.nRejectedLoops     = nnz(~keep);
    lp.fracLoopsKept      = numel(loopA) / max(numel(Apix),1);
    lp.loopDensity_mm2    = numel(loopA) / max(hullArea_mm2, eps);
    lp.meanLoopArea_um2   = safe(@mean,   loopA);
    lp.medianLoopArea_um2 = safe(@median, loopA);
    lp.sdLoopArea_um2     = safe(@std,    loopA);
    lp.iqrLoopArea_um2    = safe(@iqr,    loopA);
    lp.p95LoopArea_um2    = safe(@(x) prctile(x,95), loopA);
    lp.maxLoopArea_um2    = safe(@max,    loopA);
    lp.cvLoopArea         = safe(@std, loopA) / max(safe(@mean, loopA), eps);
    lp.totalLoopArea_mm2  = sum(loopA) / 1e6;
    lp.loopAreaFrac       = (sum(loopA)/1e6) / max(hullArea_mm2, eps);
end

function v = safe(fn, x)
    if isempty(x), v = NaN; else, v = fn(double(x)); end
end

function L = skelLength(BW)
    if ~any(BW(:)), L = 0; return; end
    right = BW(:,1:end-1) & BW(:,2:end);
    down  = BW(1:end-1,:) & BW(2:end,:);
    d1    = BW(1:end-1,1:end-1) & BW(2:end,2:end);
    d2    = BW(1:end-1,2:end)   & BW(2:end,1:end-1);
    L = nnz(right) + nnz(down) + sqrt(2)*(nnz(d1) + nnz(d2));
end