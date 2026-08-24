function T = statsFromQuRVA(vascDir, capP)
% statsFromQuRVA  Run capillaryMetrics on QuRVA's vessel segmentations.
%
%   T = statsFromQuRVA('/path/to/images/VasculatureImages')
%   T = statsFromQuRVA(vascDir, capP)
%
% QuRVA writes its vessel segmentations to a VasculatureImages folder inside
% your image folder. This reads each one, skeletonises it, and computes the
% same capillary metrics you were getting from the AutoTube masks.
%
% capP is the capillaryMetrics params struct: radiusThreshPix, minBranchPix,
% pixSzeUM, gridN. RE-TUNE radiusThreshPix for these images - whole mounts are
% usually at a different scale than your leaflet crops, so the pixel cutoff
% that separated capillaries from arterioles before will be wrong here.

    if nargin < 2, capP = struct(); end
    if ~isfield(capP,'radiusThreshPix'), capP.radiusThreshPix = 6;   end
    if ~isfield(capP,'minBranchPix'),    capP.minBranchPix    = 8;   end
    if ~isfield(capP,'pixSzeUM'),        capP.pixSzeUM        = 1.0; end
    if ~isfield(capP,'gridN'),           capP.gridN           = 8;   end

    spurLenPix = 15;    % skeleton spur pruning, via bwskel MinBranchLength

    qcDir = fullfile(fileparts(vascDir), 'capillary_qc');
    if ~exist(qcDir,'dir'), mkdir(qcDir); end

    files = [dir(fullfile(vascDir,'*.tif')); dir(fullfile(vascDir,'*.tiff')); ...
             dir(fullfile(vascDir,'*.png'))];
    if isempty(files)
        error('No images found in %s', vascDir);
    end
    fprintf('Found %d segmentations in %s\n', numel(files), vascDir);

    capTemplate = capillaryMetrics();
    capFields   = fieldnames(capTemplate);

    core = struct('fileName','', 'ok',false, 'msg','', ...
                  'maskPixels',NaN, 'skelPixels',NaN, 'maskFrac',NaN);
    for f = 1:numel(capFields), core.(capFields{f}) = NaN; end
    results = repmat(core, numel(files), 1);

    for ii = 1:numel(files)

        [~, base, ~] = fileparts(files(ii).name);
        results(ii).fileName = base;

        try
            raw = imread(fullfile(files(ii).folder, files(ii).name));

            % ---------- coerce to a binary mask ----------
            % QuRVA may write a binary mask or a colour overlay depending on
            % version and settings. Report what we got, then binarise.
            if size(raw,3) == 3
                warning(['%s is RGB - probably an overlay, not a raw mask. ' ...
                         'Binarising on any non-zero channel; check the QC image.'], base);
                im_clean = any(raw > 0, 3);
            else
                u = unique(raw(:));
                if numel(u) > 2
                    fprintf('  (%s: %d grey levels, thresholding at >0)\n', base, numel(u));
                end
                im_clean = raw > 0;
            end

            % Drop single stray pixels that would otherwise become endpoints
            im_clean = bwareaopen(im_clean, 5, 8);

            % ---------- skeleton ----------
            % bwskel is native and fast, and MinBranchLength prunes spurs in
            % the same pass. NOTE: this is NOT the same skeletonisation
            % AutoTube used, so these numbers are not directly comparable to
            % your earlier AutoTube CSV. Don't pool the two.
            im_skel = bwskel(im_clean, 'MinBranchLength', spurLenPix);

            % ---------- metrics ----------
            capPi = capP;
            capPi.hull = [];        % whole mount: convex hull is a fair proxy
            [cap, capQC] = capillaryMetrics(im_clean, im_skel, capPi);

            imwrite(capQC, fullfile(qcDir, [base '_cap.png']));

            for f = 1:numel(capFields)
                results(ii).(capFields{f}) = cap.(capFields{f});
            end
            results(ii).ok         = true;
            results(ii).maskPixels = nnz(im_clean);
            results(ii).skelPixels = nnz(im_skel);
            results(ii).maskFrac   = nnz(im_clean) / numel(im_clean);

            fprintf('[%3d/%3d] %-35s ok  junc %5d  holes %5d  meanBr %6.1f  mask %.1f%%\n', ...
                ii, numel(files), base, cap.nJunctions, cap.holeCount, ...
                cap.meanBranchLen_um, 100*results(ii).maskFrac);

        catch ME
            results(ii).ok  = false;
            results(ii).msg = ME.message;
            fprintf('[%3d/%3d] %-35s FAILED: %s\n', ii, numel(files), base, ME.message);
        end
    end

    T = struct2table(results);
    csvPath = fullfile(fileparts(vascDir), 'qurva_capillary_stats.csv');
    writetable(T, csvPath);

    fprintf('\n%d/%d succeeded.\nStats: %s\nQC: %s\n', ...
        sum([results.ok]), numel(files), csvPath, qcDir);
end