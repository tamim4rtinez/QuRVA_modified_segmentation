function exportMask(masterFolder, fileName, params)
% exportMask  Write QuRVA's vessel mask out for comparison with another
%             segmentation.
%
%   exportMask(masterFolder, 'ACI_1_cage_12_Lectin(488).tif.mat')
%   exportMask(masterFolder, fname, struct('zoomUM', 500, 'upsample', true))
%
% Writes to masterFolder/MaskExport/<base>/ :
%   mask.tif            the vessel mask, 8-bit, 0/255
%   mask_overlay.tif    mask outline drawn on the raw image
%   retina_mask.tif     the hand-drawn retinal hull
%   zoom_raw.png        small region, raw image at native resolution
%   zoom_mask.png       same region, mask
%   zoom_overlay.png    same region, outline on raw  <- the useful one
%
% THE COMPARISON THAT MATTERS is zoom_overlay at a few hundred microns
% across. At whole-retina zoom almost any segmentation looks plausible.
% What you want to see is whether the mask traces individual capillaries
% or fuses neighbouring ones into a single band.
%
% params:
%   pixSzeUM   um per pixel of the TIFF. Default 3/1.4491
%   zoomUM     width of the zoom crop. Default 400
%   zoomAt     [x y] centre of the crop in MASK pixel coords. Default: a
%              dense region found automatically
%   upsample   true writes the mask at the original TIFF resolution rather
%              than QuRVA's downsampled one, so it can be overlaid directly
%              on another pipeline's mask. Default false.

    if nargin < 3, params = struct(); end
    p.pixSzeUM = getf(params,'pixSzeUM', 3/1.4491);
    p.zoomUM   = getf(params,'zoomUM',   400);
    p.zoomAt   = getf(params,'zoomAt',   []);
    p.upsample = getf(params,'upsample', false);

    S = load(fullfile(masterFolder,'VasculatureNumbers',fileName), 'smoothVessels');
    M = load(fullfile(masterFolder,'Masks',fileName), 'thisMask');

    im_clean = logical(S.smoothVessels);
    im_hull  = logical(M.thisMask);

        sf = 1;
    if ~isequal(size(im_hull), size(im_clean))
        sf = size(im_hull,1)/size(im_clean,1);
        hull_full = im_hull;                              % keep for output
        im_hull   = imresize(im_hull, size(im_clean), 'nearest');
    else
        hull_full = im_hull;
    end
    pxUM = p.pixSzeUM * sf;            % um per pixel of im_clean

    base   = erase(fileName, '.mat');
    outDir = fullfile(masterFolder, 'MaskExport', matlab.lang.makeValidName(base));
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % ---- raw image, if it is still beside the .mat files ----
    raw = [];
    rawPath = fullfile(masterFolder, base);
    if exist(rawPath, 'file')
        R = imread(rawPath);
        raw = im2uint8(R(:,:,1));                    % QuRVA uses the red channel
        raw = imresize(raw, size(im_clean));         % match mask resolution
    else
        fprintf('Raw image not found at %s - overlays will be skipped.\n', rawPath);
    end

    % ---- whole-retina exports ----
    if p.upsample
        outMask = imresize(im_clean, size(hull_full), 'nearest');
        fprintf('Writing mask upsampled to %dx%d (original TIFF size).\n', ...
                size(outMask,1), size(outMask,2));
    else
        outMask = im_clean;
    end
    imwrite(uint8(outMask)*255, fullfile(outDir,'mask.tif'));
        imwrite(uint8(hull_full)*255, fullfile(outDir,'retina_mask.tif'));
    if ~isempty(raw)
        ov = cat(3, raw, raw, raw);
        edge = bwperim(im_clean);
        ovR = ov(:,:,1); ovG = ov(:,:,2); ovB = ov(:,:,3);
        ovR(edge) = 0; ovG(edge) = 255; ovB(edge) = 0;
        ov = cat(3, ovR, ovG, ovB);
        imwrite(ov, fullfile(outDir,'mask_overlay.tif'));
    end

    % ---- pick a zoom region ----
    z = max(round(p.zoomUM / pxUM), 32);
    if isempty(p.zoomAt)
        % densest z-by-z window inside the retina
        dens = conv2(double(im_clean & im_hull), ones(z)/z^2, 'same');
        dens(~im_hull) = 0;
        [~, idx] = max(dens(:));
        [cy, cx] = ind2sub(size(dens), idx);
    else
        cx = p.zoomAt(1); cy = p.zoomAt(2);
    end
    r0 = max(min(round(cy - z/2), size(im_clean,1)-z), 1);
    c0 = max(min(round(cx - z/2), size(im_clean,2)-z), 1);
    rr = r0:r0+z-1; cc = c0:c0+z-1;

    imwrite(uint8(im_clean(rr,cc))*255, fullfile(outDir,'zoom_mask.png'));
    if ~isempty(raw)
        imwrite(raw(rr,cc), fullfile(outDir,'zoom_raw.png'));
        sub  = raw(rr,cc);
        edge = bwperim(im_clean(rr,cc));
        zR = sub; zG = sub; zB = sub;
        zR(edge) = 0; zG(edge) = 255; zB(edge) = 0;
        imwrite(cat(3,zR,zG,zB), fullfile(outDir,'zoom_overlay.png'));
    end

    % ---- numbers to quote alongside the images ----
    D = bwdist(~im_clean);
    sk = bwmorph(im_clean, 'thin', Inf);
    fprintf('\n%s\n', base);
    fprintf('  mask resolution      %d x %d\n', size(im_clean,1), size(im_clean,2));
    fprintf('  TIFF resolution      %d x %d\n', size(im_hull,1),  size(im_hull,2));
    fprintf('  scaleFactor          %.3f\n', sf);
    fprintf('  um per mask pixel    %.3f\n', pxUM);
    fprintf('  vessel area fraction %.4f  (vessel px / retina px)\n', ...
            nnz(im_clean & im_hull) / max(nnz(im_hull),1));
    fprintf('  mean vessel radius   %.2f um   <- real capillaries are 2-3 um\n', ...
            mean(D(sk)) * pxUM);
    fprintf('  zoom centred at      [%d %d], %.0f um across\n', cx, cy, z*pxUM);
    fprintf('\nWrote %s\n\n', outDir);
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end