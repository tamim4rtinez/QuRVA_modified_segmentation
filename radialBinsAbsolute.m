function T = radialBinsAbsolute(masterFolder, params)
% radialBinsAbsolute  Capillary density in FIXED-WIDTH annuli from the ONH.
%
%   T = radialBinsAbsolute(masterFolder)
%   T = radialBinsAbsolute(masterFolder, params)
%
% Written to match a manual workflow that measures a fixed distance from the
% optic nerve head (e.g. 0-250, 250-500, 500-750 um) rather than scaling the
% bins to each retina. Use this ONLY for the cross-validation against manual
% counts; the main pipeline uses relative (normalised-eccentricity) bins.
%
% params fields (all optional):
%   basePxUM         um per pixel of the TIFF that was segmented.
%                    Series-2 mouse: 3/2.899 = 1.0348.   Default 1.0348.
%   binEdges_um      annulus edges from the ONH. Default [0 250 500 750].
%   radiusThreshUM   capillary/large-vessel cutoff, RADIUS in um. Default 5.
%   minBranchPix     spur removal, in pixels. Default 8.
%   saveQC           write a PNG per retina showing the rings. Default true.
%
% Output columns, one row per retina per bin:
%   fileName, bin, inner_um, outer_um,
%   tissueArea_mm2      area of retina inside that annulus
%   capLen_mm           capillary skeleton length in that annulus
%   lenDensity_mm_mm2   capLen_mm / tissueArea_mm2      <- length density
%   areaFrac            capillary MASK pixels / tissue pixels  <- area fraction
%   coverage            fraction of the annulus that is inside the retina
%                       (< 1 means the ring runs off the tissue edge)
%
% Both lenDensity and areaFrac are reported because "capillary density" is
% measured either way in the literature. Compare against whichever your
% colleague computed - they are NOT interchangeable.

    if nargin < 2, params = struct(); end
    p.basePxUM       = getf(params, 'basePxUM',       3/2.899);
    p.binEdges_um    = getf(params, 'binEdges_um',    [0 250 500 750]);
    p.radiusThreshUM = getf(params, 'radiusThreshUM', 5);
    p.minBranchPix   = getf(params, 'minBranchPix',   8);
    p.saveQC         = getf(params, 'saveQC',         true);

    vascDir = fullfile(masterFolder, 'VasculatureNumbers');
    maskDir = fullfile(masterFolder, 'Masks');
    ocDir   = fullfile(masterFolder, 'ONCenter');
    qcDir   = fullfile(masterFolder, 'RadialQC');

    assert(exist(vascDir,'dir') == 7, ...
        'No VasculatureNumbers folder in %s', masterFolder);
    if p.saveQC && ~exist(qcDir,'dir'), mkdir(qcDir); end

    files  = dir(fullfile(vascDir, '*.mat'));
    nBins  = numel(p.binEdges_um) - 1;
    assert(~isempty(files), 'No .mat files in %s', vascDir);
    fprintf('Found %d retinas, %d bins: %s um\n', ...
        numel(files), nBins, mat2str(p.binEdges_um));

    rows = struct('fileName',{},'bin',{},'inner_um',{},'outer_um',{}, ...
                  'tissueArea_mm2',{},'capLen_mm',{}, ...
                  'lenDensity_mm_mm2',{},'areaFrac',{},'coverage',{});

    for ii = 1:numel(files)

        base = erase(files(ii).name, '.mat');

        try
            S = load(fullfile(vascDir, files(ii).name), ...
                     'smoothVessels', 'vesselSkelMask');
            M = load(fullfile(maskDir, files(ii).name), 'thisMask');
            C = load(fullfile(ocDir,   files(ii).name), 'thisONCenter');

            im_clean = logical(S.smoothVessels);
            im_skel  = logical(S.vesselSkelMask);
            im_hull  = logical(M.thisMask);

            % QuRVA rescales before segmenting; match mask to vessel space
            scaleFactor = 1;
            if ~isequal(size(im_hull), size(im_clean))
                scaleFactor = size(im_hull,1) / size(im_clean,1);
                im_hull = imresize(im_hull, size(im_clean), 'nearest');
            end
            apexXY = C.thisONCenter(:).' / scaleFactor;   % [x y]

            pxUM   = p.basePxUM * scaleFactor;   % um per pixel of im_clean
            pxArea = pxUM^2;

            % ---- capillary / large-vessel split, threshold in microns ----
            threshPix = p.radiusThreshUM / pxUM;
            D = bwdist(~im_clean);
            skelRadius = zeros(size(D));
            skelRadius(im_skel) = D(im_skel);
            skel_cap = im_skel & skelRadius > 0 & skelRadius <= threshPix;
            skel_cap = bwareaopen(skel_cap, p.minBranchPix, 8);
            mask_cap = im_clean & (D <= threshPix);

            % ---- distance from the ONH, in microns ----
            [nr, nc] = size(im_hull);
            [X, Y]   = meshgrid(1:nc, 1:nr);
            Dap_um   = hypot(X - apexXY(1), Y - apexXY(2)) * pxUM;

            for b = 1:nBins
                lo = p.binEdges_um(b);
                hi = p.binEdges_um(b+1);

                ring   = Dap_um >= lo & Dap_um < hi;   % full annulus
                band   = im_hull & ring;               % annulus inside retina

                tissuePix = nnz(band);
                ringPix   = nnz(ring);

                skelPix = nnz(skel_cap & band);
                maskPix = nnz(mask_cap & band);
                r.fileName          = base;
                r.bin               = b;
                r.inner_um          = lo;
                r.outer_um          = hi;
                r.tissueArea_mm2    = tissuePix * pxArea / 1e6;
                r.capLen_mm         = skelPix   * pxUM   / 1e3;
                r.lenDensity_mm_mm2 = r.capLen_mm / max(r.tissueArea_mm2, eps);
                r.areaFrac          = maskPix / max(tissuePix, 1);
                r.coverage          = tissuePix / max(ringPix, 1);

                rows(end+1) = r; %#ok<AGROW>
            end

            fprintf(['[%3d/%3d] %-40s  px %.4f um  ' ...
                     'dens:'], ii, numel(files), base, pxUM);
            fprintf(' %6.2f', [rows(end-nBins+1:end).lenDensity_mm_mm2]);
            fprintf('   cov:');
            fprintf(' %4.2f', [rows(end-nBins+1:end).coverage]);
            fprintf('\n');

            % ---- QC overlay ----
            if p.saveQC
                rgb = repmat(0.18*double(im_clean), [1 1 3]);
                rgb(:,:,2) = max(rgb(:,:,2), ...
                    double(imdilate(skel_cap, strel('disk',1))));
                                    skel_large = im_skel & skelRadius > threshPix;
                lrg = double(imdilate(skel_large, strel('disk',1)));
                rgb(:,:,1) = max(rgb(:,:,1), lrg);
                rgb(:,:,3) = max(rgb(:,:,3), lrg);
                f = figure('Visible','off','Color','w');
                imshow(rgb); hold on
                th = linspace(0, 2*pi, 400);
                for b = 2:numel(p.binEdges_um)
                    rp = p.binEdges_um(b) / pxUM;
                    plot(apexXY(1)+rp*cos(th), apexXY(2)+rp*sin(th), ...
                         'y-', 'LineWidth', 1.2);
                end
                plot(apexXY(1), apexXY(2), 'y*', 'MarkerSize', 12);
                title(strrep(base,'_','\_'), 'Interpreter','tex');
                exportgraphics(f, fullfile(qcDir, [base '_rings.png']), ...
                               'Resolution', 120);
                close(f);
            end

        catch ME
            fprintf('[%3d/%3d] %-40s  FAILED: %s\n', ...
                ii, numel(files), base, ME.message);
        end
    end

    T = struct2table(rows);
    outPath = fullfile(masterFolder, 'radial_bins_absolute.csv');
    writetable(T, outPath);
    fprintf('\nWrote %s\nQC rings: %s\n', outPath, qcDir);
end

% ----------------------------------------------------------------------
function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end