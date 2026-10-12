function exportSkeleton(masterFolder, fileName, params)
% exportSkeleton  Write QuRVA's vessel skeleton out for overlay in Fiji.
%
%   exportSkeleton(masterFolder, 'ACI_1_cage_12_Lectin(488).tif.mat')
%   exportSkeleton(masterFolder, fname, struct('upsample', true, 'thickenPix', 2))
%
% Writes to masterFolder/MaskExport/<base>/ :
%   skeleton.tif          full skeleton, all vessels
%   skeleton_cap.tif      capillary skeleton only (radius <= radiusThreshUM)
%   skeleton_large.tif    large-vessel skeleton only
%   skeleton_rgb.tif      green = capillary, magenta = large, on black
%   zoom_skel_overlay.png skeleton on the raw image, small region
%
% UPSAMPLING: a one-pixel skeleton resampled to a larger grid breaks into
% dashes, because nearest-neighbour cannot interpolate a line. So when
% upsample is on, the skeleton is dilated by thickenPix FIRST and then
% resized. The result is a continuous line you can see at full zoom, but it
% is a display object - do not measure length from the upsampled file.
% Measure from skeleton.tif at native resolution, or from MATLAB.
%
% params:
%   pixSzeUM        um per pixel of the TIFF. Default 3/1.4491
%   radiusThreshUM  capillary / large cutoff, radius in um. Default 17.1
%   minBranchPix    spur removal for the capillary skeleton. Default 8
%   upsample        write at the original TIFF resolution. Default false
%   thickenPix      dilation for visibility. Default 1 (2-3 reads better
%                   when overlaid on a large image)
%   zoomUM          width of the zoom crop. Default 400
%   zoomAt          [x y] crop centre in MASK pixel coords. Default: auto

    if nargin < 3, params = struct(); end
    p.pixSzeUM       = getf(params,'pixSzeUM',       3/1.4491);
    p.radiusThreshUM = getf(params,'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params,'minBranchPix',   8);
    p.upsample       = getf(params,'upsample',       false);
    p.thickenPix     = getf(params,'thickenPix',     1);
    p.zoomUM         = getf(params,'zoomUM',         400);
    p.zoomAt         = getf(params,'zoomAt',         []);

    S = load(fullfile(masterFolder,'VasculatureNumbers',fileName), ...
             'smoothVessels','vesselSkelMask');
    M = load(fullfile(masterFolder,'Masks',fileName), 'thisMask');

    im_clean  = logical(S.smoothVessels);
    im_skel   = logical(S.vesselSkelMask);
    hull_full = logical(M.thisMask);

    sf = 1;
    if ~isequal(size(hull_full), size(im_clean))
        sf = size(hull_full,1)/size(im_clean,1);
    end
    pxUM = p.pixSzeUM * sf;

    % ---- capillary / large split ----
    D  = bwdist(~im_clean);
    sr = zeros(size(D)); sr(im_skel) = D(im_skel);
    threshPix  = p.radiusThreshUM / pxUM;
    skel_cap   = bwareaopen(im_skel & sr > 0 & sr <= threshPix, p.minBranchPix, 8);
    skel_large = im_skel & sr > threshPix;

    base   = erase(fileName, '.mat');
    outDir = fullfile(masterFolder,'MaskExport', matlab.lang.makeValidName(base));
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % ---- write, thickened and optionally upsampled ----
    tgt = size(hull_full);
    write_one(im_skel,    'skeleton.tif');
    write_one(skel_cap,   'skeleton_cap.tif');
    write_one(skel_large, 'skeleton_large.tif');

    % colour composite: green capillary, magenta large
    g = prep(skel_cap); m = prep(skel_large);
    rgb = cat(3, uint8(m)*255, uint8(g)*255, uint8(m)*255);
    imwrite(rgb, fullfile(outDir,'skeleton_rgb.tif'));

    % ---- zoom overlay on the raw image ----
    raw = [];
    rawPath = fullfile(masterFolder, base);
    if exist(rawPath,'file')
        R = imread(rawPath);
        raw = imresize(im2uint8(R(:,:,1)), size(im_clean));
    end

    z = max(round(p.zoomUM / pxUM), 32);
    if isempty(p.zoomAt)
        dens = conv2(double(im_clean), ones(z)/z^2, 'same');
        [~, idx] = max(dens(:));
        [cy, cx] = ind2sub(size(dens), idx);
    else
        cx = p.zoomAt(1); cy = p.zoomAt(2);
    end
    r0 = max(min(round(cy-z/2), size(im_clean,1)-z), 1);
    c0 = max(min(round(cx-z/2), size(im_clean,2)-z), 1);
    rr = r0:r0+z-1; cc = c0:c0+z-1;

    if ~isempty(raw)
        sub = raw(rr,cc);
        sc  = skel_cap(rr,cc);  sl = skel_large(rr,cc);
        zR = sub; zG = sub; zB = sub;
        zR(sc) = 0;   zG(sc) = 255; zB(sc) = 0;
        zR(sl) = 255; zG(sl) = 0;   zB(sl) = 255;
        imwrite(cat(3,zR,zG,zB), fullfile(outDir,'zoom_skel_overlay.png'));
    end

    fprintf('\n%s\n', base);
    fprintf('  um per mask pixel  %.3f     scaleFactor %.3f\n', pxUM, sf);
    fprintf('  skeleton pixels    %d total, %d capillary, %d large\n', ...
            nnz(im_skel), nnz(skel_cap), nnz(skel_large));
    fprintf('  capillary cutoff   %.1f um radius (%.2f px)\n', ...
            p.radiusThreshUM, threshPix);
    if p.upsample
        fprintf('  written upsampled to %d x %d, dilated %d px for visibility\n', ...
                tgt(1), tgt(2), p.thickenPix);
        fprintf('  -> DISPLAY ONLY. Do not measure length from these files.\n');
    end
    fprintf('\nWrote %s\n\n', outDir);

    % ------------------------------------------------------------------
    function out = prep(BW)
        if p.thickenPix > 0
            BW = imdilate(BW, strel('disk', p.thickenPix));
        end
        if p.upsample
            BW = imresize(BW, tgt, 'nearest');
        end
        out = BW;
    end

    function write_one(BW, name)
        imwrite(uint8(prep(BW))*255, fullfile(outDir, name));
    end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end