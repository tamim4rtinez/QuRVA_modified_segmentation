function exportRadialROIs(masterFolder, fileName, params)
% exportRadialROIs  Write per-bin masks for OrientationJ analysis in FIJI.
%
%   exportRadialROIs(masterFolder, 'ACI_1_cage_12_Lectin(488).tif.mat')
%
% Uses the SAME radial binning as radialTopology: equal-width bins spanning
% each retina's own min-to-max distance from the optic nerve head.
%
% Writes to masterFolder/OrientationJ/<base>/ :
%   capillary_mask.tif   8-bit, capillaries only (large vessels removed)
%   band_1.tif ... band_N.tif   8-bit band masks, already clipped to retina
%   bin_edges.txt        bin radii in um, for the record
%
% In FIJI: open band_k.tif -> Edit > Selection > Create Selection,
% then select capillary_mask.tif -> Edit > Selection > Restore Selection,
% then Analyze > OrientationJ > OrientationJ Measure.

    if nargin < 3, params = struct(); end
    p.basePxUM       = getf(params, 'basePxUM',       3/1.4491);
    p.nBins          = getf(params, 'nBins',          6);
    p.radiusThreshUM = getf(params, 'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params, 'minBranchPix',   8);

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

    % ---- capillary-only MASK (not skeleton - OrientationJ needs width) ----
    D = bwdist(~im_clean);
    skelRadius = zeros(size(D));
    skelRadius(im_skel) = D(im_skel);
    threshPix = p.radiusThreshUM / pxUM;
    skel_cap  = bwareaopen(im_skel & skelRadius > 0 & skelRadius <= threshPix, ...
                           p.minBranchPix, 8);
    mask_cap  = im_clean & (D <= threshPix);
    mask_cap  = mask_cap & imdilate(skel_cap, strel('disk', ceil(threshPix)));

    % ---- radial bands, identical to radialTopology ----
    [nr, nc] = size(im_hull);
    [X, Y]   = meshgrid(1:nc, 1:nr);
    Dap      = hypot(X - apexXY(1), Y - apexXY(2));
    dHull    = Dap(im_hull);
    edges    = linspace(min(dHull), max(dHull)+eps, p.nBins+1);

    base   = erase(fileName, '.mat');
    outDir = fullfile(masterFolder, 'OrientationJ', matlab.lang.makeValidName(base));
    if ~exist(outDir,'dir'), mkdir(outDir); end

    imwrite(uint8(mask_cap)*255, fullfile(outDir,'capillary_mask.tif'));

    fid = fopen(fullfile(outDir,'bin_edges.txt'),'w');
    fprintf(fid, 'file: %s\npixel size: %.4f um\napex: %.1f %.1f\n\n', ...
            base, pxUM, apexXY(1), apexXY(2));
    fprintf('\n%s   px %.4f um\n', base, pxUM);

    for b = 1:p.nBins
        band = im_hull & Dap >= edges(b) & Dap < edges(b+1);
        imwrite(uint8(band)*255, fullfile(outDir, sprintf('band_%d.tif', b)));

        lo = edges(b)*pxUM;  hi = edges(b+1)*pxUM;
        vf = nnz(mask_cap & band) / max(nnz(band),1);
        fprintf(fid, 'band %d: %8.1f - %8.1f um   vessel fraction %.3f\n', b, lo, hi, vf);
        fprintf('  band %d: %7.0f - %7.0f um   vessel fraction %.3f\n', b, lo, hi, vf);
    end
    fclose(fid);

    fprintf('\nWrote %s\n', outDir);
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end