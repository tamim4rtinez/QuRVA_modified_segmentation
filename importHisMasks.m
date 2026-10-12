function T = importHisMasks(hisDir, oldMaster, newMaster, params)
% importHisMasks  Build a QuRVA-style analysis folder from externally
%                 produced vessel masks and retinal hulls.
%
%   T = importHisMasks(hisDir, oldMaster, newMaster)
%
% Pairs, by strain and animal number:
%   *_whole_vessel_mask_*.tif    -> smoothVessels   (the vessel mask)
%   *_analysis_mask_75um_*.tif   -> thisMask        (the retinal hull)
% Both come from the same pipeline in the same coordinate space, so no
% registration is needed.
%
% The skeleton is made the way QuRVA makes it - bwmorph 'thin' to
% convergence - so skeletonisation is not a variable between this analysis
% and your original one. QuRVA's two extra steps are deliberately skipped:
%
%   - the pre-skeleton imdilate(..., DilatingRadius). That dilation is why
%     the QuRVA masks are ~1.5x wider than these and why neighbouring
%     capillaries fuse. These masks are used as given.
%   - the retinal-edge strip (imdilate(bwperim(hull), disk 5)). The 75um
%     hull is already inset from the tissue edge.
%
% Output names follow your existing convention, Strain_N_cage_M, with the
% cage number taken from the matching QuRVA filename. Where no cage number
% can be found the name falls back to Strain_N and that is flagged in the log.
%
% Writes:
%   newMaster/VasculatureNumbers/<name>.mat   smoothVessels, vesselSkelMask,
%                                             brchPts, endPts, aVascZone
%   newMaster/Masks/<name>.mat                thisMask
%   newMaster/import_log.csv
%
% ONCenter/ is NOT written here - run pickONH(newMaster) afterwards.
%
% ---------------------------------------------------------------------
% RE-TUNE BEFORE ANALYSING. Both of these were set for dilated QuRVA masks
% at 3.79 um/pixel and will be wrong here:
%
%   radiusThreshUM  the capillary / large-vessel cutoff. 17.1 um was derived
%                   from a pixel threshold on dilated masks; on an
%                   undilated segmentation it will classify almost
%                   everything as capillary. The log prints skeleton radius
%                   percentiles - set it from those, around p90-p95.
%   minBranchPix    spur pruning in PIXELS. Convert: <um> / pixSzeUM.
% ---------------------------------------------------------------------
%
% params:
%   pixSzeUM    um per pixel of his masks at full size. Default 1/1.4491
%   workScale   fraction of full resolution to work at. Default 1.
%               bwmorph 'thin' on a 400+ Mpx image is slow and memory
%               heavy; 0.5 is much faster and still ~2.7x finer than the
%               QuRVA masks.
%   overwrite   redo retinas already imported. Default false

    if nargin < 4, params = struct(); end
    p.pixSzeUM  = getf(params,'pixSzeUM',  1/1.4491);
    p.workScale = getf(params,'workScale', 1);
    p.overwrite = getf(params,'overwrite', false);

    pxUM = p.pixSzeUM / p.workScale;

    for d = {'VasculatureNumbers','Masks'}
        dd = fullfile(newMaster, d{1});
        if ~exist(dd,'dir'), mkdir(dd); end
    end

    vessFiles = dir(fullfile(hisDir,'*whole_vessel_mask*.tif'));
    hullFiles = dir(fullfile(hisDir,'*analysis_mask*.tif'));
    assert(~isempty(vessFiles), 'No *whole_vessel_mask*.tif in %s', hisDir);
    assert(~isempty(hullFiles), 'No *analysis_mask*.tif in %s', hisDir);

    keyOf = @(n) upper(char(join(regexp(n,'^([A-Za-z0-9]+)_(\d+)_','tokens','once'),'_')));
    vKey = cellfun(@(n) keyOf(n), {vessFiles.name}, 'UniformOutput', false);
    hKey = cellfun(@(n) keyOf(n), {hullFiles.name}, 'UniformOutput', false);

    % cage numbers from your existing QuRVA filenames
    oldF = dir(fullfile(oldMaster,'VasculatureNumbers','*.mat'));
    cageOf = containers.Map('KeyType','char','ValueType','char');
    for k = 1:numel(oldF)
        t = regexp(oldF(k).name,'^([A-Za-z0-9]+)_(\d+)','tokens','once');
        if isempty(t), continue; end
        c = regexp(oldF(k).name,'cage_?(\d+)','tokens','once');
        kk = upper([t{1} '_' t{2}]);
        if ~isempty(c), cageOf(kk) = c{1}; end
    end

    fprintf('\nImporting at %.3f um/pixel (%.0f%% of full resolution)\n\n', ...
            pxUM, 100*p.workScale);

    rows = [];
    for ii = 1:numel(vessFiles)
        key = vKey{ii};
        if isempty(key)
            fprintf('  SKIP  %s - cannot parse strain/animal\n', vessFiles(ii).name);
            continue
        end
        jj = find(strcmp(hKey, key), 1);
        if isempty(jj)
            fprintf('  SKIP  %s - no analysis mask for %s\n', vessFiles(ii).name, key);
            continue
        end

        if isKey(cageOf, key)
            outBase = sprintf('%s_cage_%s', key, cageOf(key));
            cageNote = "";
        else
            outBase = key;
            cageNote = "no cage number found";
        end
        outName = [outBase '.mat'];

        if ~p.overwrite && exist(fullfile(newMaster,'VasculatureNumbers',outName),'file')
            fprintf('  have  %s\n', outBase); continue
        end

        try
            V = readBin(fullfile(vessFiles(ii).folder, vessFiles(ii).name), p.workScale);
            H = readBin(fullfile(hullFiles(jj).folder, hullFiles(jj).name), p.workScale);
            if ~isequal(size(V), size(H))
                H = imresize(H, size(V), 'nearest');
            end

            smoothVessels = V & H;

            % QuRVA's skeleton step, without the edge strip
            vesselSkelMask = bwmorph(smoothVessels, 'thin', Inf) & smoothVessels;
  spurPix = round(getf(params,'spurUM',20) / pxUM);
            brch = bwmorph(vesselSkelMask,'branchpoints');
chopped = vesselSkelMask & ~brch;
chopped = bwareaopen(chopped, spurPix);
vesselSkelMask = bwareaopen(chopped | brch, 1);

            brchPts   = bwmorph(vesselSkelMask, 'branchpoints');
            endPts    = bwmorph(vesselSkelMask, 'endpoints');
            aVascZone = false(size(smoothVessels));
            thisMask  = H;                                       %#ok<NASGU>

            D  = bwdist(~smoothVessels);
            rr = D(vesselSkelMask) * pxUM;
            q  = prctile(rr, [50 75 90 95 99]);

            save(fullfile(newMaster,'VasculatureNumbers',outName), ...
                 'smoothVessels','vesselSkelMask','brchPts','endPts', ...
                 'aVascZone','-v7.3');
            save(fullfile(newMaster,'Masks',outName), 'thisMask','-v7.3');

            r.key = string(key); r.strain = string(extractBefore(key,'_'));
            r.outName = string(outBase); r.note = cageNote;
            r.vesselFile = string(vessFiles(ii).name);
            r.hullFile   = string(hullFiles(jj).name);
            r.height = size(V,1); r.width = size(V,2); r.pxSzeUM = pxUM;
            r.hullArea_mm2 = nnz(H) * pxUM^2 / 1e6;
            r.vesselFrac   = nnz(smoothVessels) / max(nnz(H),1);
            r.skelPix      = nnz(vesselSkelMask);
            r.radius_p50 = q(1); r.radius_p75 = q(2);
            r.radius_p90 = q(3); r.radius_p95 = q(4); r.radius_p99 = q(5);
            rows = [rows; r]; %#ok<AGROW>

            fprintf(['  %-22s  %.1f mm2  vesselFrac %.3f  ' ...
                     'radius um p50 %.2f p90 %.2f p99 %.2f  %s\n'], ...
                outBase, r.hullArea_mm2, r.vesselFrac, q(1), q(3), q(5), cageNote);

            clear V H smoothVessels vesselSkelMask brchPts endPts aVascZone thisMask D

        catch ME
            fprintf('  FAIL  %s: %s\n', key, ME.message);
        end
    end

    assert(~isempty(rows), 'Nothing imported.');
    T = struct2table(rows);
    writetable(T, fullfile(newMaster,'import_log.csv'));

    fprintf('\nImported %d retinas to %s\n', height(T), newMaster);
    fprintf('\nNEXT: pickONH(''%s'')  to set the optic nerve head positions.\n', newMaster);
    fprintf('\nTHEN set:\n');
    fprintf('  capP.pixSzeUM       = %.6f;\n', pxUM);
    fprintf('  capP.radiusThreshUM = <around p90-p95 below>;\n');
    fprintf('  capP.minBranchPix   = round(<um> / %.4f);\n\n', pxUM);
    fprintf('Median across retinas:  p50 %.2f   p90 %.2f   p95 %.2f   p99 %.2f um\n\n', ...
        median(T.radius_p50), median(T.radius_p90), ...
        median(T.radius_p95), median(T.radius_p99));
end

% ======================================================================
function BW = readBin(path, scale)
    I = imread(path);
    if size(I,3) > 1, BW = any(I > 0, 3); else, BW = I > 0; end
    if scale < 1, BW = imresize(BW, scale, 'nearest'); end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end