function T = vesselLengthDensity(masterFolder, params)
% vesselLengthDensity  Whole-retina length density, all vessels and split.
%
%   T = vesselLengthDensity(masterFolder, struct('pixSzeUM', 3/1.4491))
%
% Written to compare directly against pipelines that do NOT exclude large
% vessels. Reports three densities from the same skeleton:
%
%   allLenDensity_mm_mm2    every vessel  <- compare this to other methods
%   capLenDensity_mm_mm2    capillaries only
%   largeLenDensity_mm_mm2  large vessels only
%
% Lengths are true path lengths along the skeleton (diagonal steps counted
% as sqrt(2)), not pixel counts - a pixel count undercounts diagonal runs
% by up to 41%.
%
% ALSO REPORTED, because these are what tell you whether a disagreement
% with another pipeline is about vessel classification or about how much
% vessel each pipeline resolves at all:
%
%   largeLenFrac      fraction of total length classified as large vessel.
%                     If this is ~0.1, excluding large vessels cannot
%                     explain more than a ~10% difference in density.
%   meanCapRadius_um  mean vessel radius along the capillary skeleton.
%                     Should be ~2-4 um for real capillaries. Much larger
%                     means the mask is merging adjacent capillaries into
%                     single bands, and total length is being undercounted
%                     by roughly the same factor.
%   skelPixelsPerMm2  skeleton pixels per mm2 of retina, no unit conversion
%                     applied - immune to pixSzeUM errors.
%
% params:
%   pixSzeUM         REQUIRED. um per pixel of the TIFF.
%   radiusThreshUM   capillary / large cutoff, radius in um. Default 17.1
%   minBranchPix     spur removal. Default 8

    if nargin < 2, params = struct(); end
    if ~isfield(params,'pixSzeUM') || isempty(params.pixSzeUM)
        error('params.pixSzeUM is required (microns per pixel of the TIFF).');
    end
    p.pixSzeUM       = params.pixSzeUM;
    p.radiusThreshUM = getf(params,'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params,'minBranchPix',   8);

    vascDir = fullfile(masterFolder,'VasculatureNumbers');
    maskDir = fullfile(masterFolder,'Masks');
    files   = dir(fullfile(vascDir,'*.mat'));
    assert(~isempty(files), 'No .mat files in %s', vascDir);

    rows = [];
    fprintf('\n%-42s %7s %7s %7s %7s %7s\n', 'file', ...
            'all', 'cap', 'large', 'lrgFrac', 'capRad');
    fprintf('%s\n', repmat('-', 1, 90));

    for ii = 1:numel(files)
        base = erase(files(ii).name, '.mat');
        try
            S = load(fullfile(vascDir, files(ii).name), ...
                     'smoothVessels','vesselSkelMask');
            M = load(fullfile(maskDir, files(ii).name), 'thisMask');

            im_clean = logical(S.smoothVessels);
            im_skel  = logical(S.vesselSkelMask);
            im_hull  = logical(M.thisMask);

            sf = 1;
            if ~isequal(size(im_hull), size(im_clean))
                sf = size(im_hull,1)/size(im_clean,1);
                im_hull = imresize(im_hull, size(im_clean), 'nearest');
            end
            pxUM = p.pixSzeUM * sf;

            D  = bwdist(~im_clean);
            sr = zeros(size(D)); sr(im_skel) = D(im_skel);
            threshPix = p.radiusThreshUM / pxUM;

            skel_cap   = bwareaopen(im_skel & sr > 0 & sr <= threshPix, ...
                                    p.minBranchPix, 8);
            skel_large = im_skel & sr > threshPix;

            area_mm2 = nnz(im_hull) * pxUM^2 / 1e6;

            allLen   = skelLength(im_skel)    * pxUM / 1e3;   % mm
            capLen   = skelLength(skel_cap)   * pxUM / 1e3;
            largeLen = skelLength(skel_large) * pxUM / 1e3;

            capRad = mean(D(skel_cap)) * pxUM;

            r.fileName   = string(base);
            r.strain     = string(extractBefore(char(base),'_'));
            r.pxUM       = pxUM;
            r.area_mm2   = area_mm2;

            r.allLen_mm   = allLen;
            r.capLen_mm   = capLen;
            r.largeLen_mm = largeLen;

            r.allLenDensity_mm_mm2   = allLen   / max(area_mm2, eps);
            r.capLenDensity_mm_mm2   = capLen   / max(area_mm2, eps);
            r.largeLenDensity_mm_mm2 = largeLen / max(area_mm2, eps);

            r.largeLenFrac     = largeLen / max(allLen, eps);
            r.meanCapRadius_um = capRad;
            r.skelPixelsPerMm2 = nnz(im_skel) / max(area_mm2, eps);
            r.vesselAreaFrac   = nnz(im_clean) / max(nnz(im_hull), 1);

            rows = [rows; r]; %#ok<AGROW>

            fprintf('%-42s %7.1f %7.1f %7.2f %7.3f %7.1f\n', base, ...
                r.allLenDensity_mm_mm2, r.capLenDensity_mm_mm2, ...
                r.largeLenDensity_mm_mm2, r.largeLenFrac, capRad);

        catch ME
            fprintf('%-42s FAILED: %s\n', base, ME.message);
        end
    end

    T = struct2table(rows);
    out = fullfile(masterFolder,'vessel_length_density.csv');
    writetable(T, out);
    fprintf('\nWrote %s\n', out);

    fprintf('\n--- by strain ---\n');
    disp(groupsummary(T, 'strain', 'mean', ...
        {'allLenDensity_mm_mm2','capLenDensity_mm_mm2', ...
         'largeLenFrac','meanCapRadius_um','vesselAreaFrac'}));

    fprintf(['\nIf largeLenFrac is small (~0.1), excluding large vessels\n' ...
             'cannot account for more than about a 10%% difference against\n' ...
             'a pipeline that includes them. A larger gap than that points\n' ...
             'at segmentation or resolution instead - check meanCapRadius_um,\n' ...
             'which should be roughly 2-4 um for resolved capillaries.\n\n']);
end

% ======================================================================
function L = skelLength(BW)
% True path length of a skeleton in pixels: orthogonal steps 1,
% diagonal steps sqrt(2). Each adjacency counted once.
    if ~any(BW(:)), L = 0; return; end

    % count adjacencies directly
    right = BW(:,1:end-1) & BW(:,2:end);
    down  = BW(1:end-1,:) & BW(2:end,:);
    dia1  = BW(1:end-1,1:end-1) & BW(2:end,2:end);
    dia2  = BW(1:end-1,2:end)   & BW(2:end,1:end-1);
    L = nnz(right) + nnz(down) + sqrt(2)*(nnz(dia1) + nnz(dia2));
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end