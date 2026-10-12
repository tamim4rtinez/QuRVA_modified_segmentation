function T = compareOneMask(masterFolder, hisDir, key, params)
% compareOneMask  Register and compare ONE retina, appending to a shared
%                 output folder. Run it again with a different key and the
%                 results accumulate.
%
%   compareOneMask(masterFolder, hisDir, 'ACI_1')
%   compareOneMask(masterFolder, hisDir, 'BN_3', struct('workScale',1))
%
% Builds up, across runs:
%   MaskComparison/overlays/<KEY>_overlay.png
%   MaskComparison/crops/<KEY>_crop1..4.png
%   MaskComparison/my_masks/<KEY>_qurva_mask_registered.tif
%   MaskComparison/mask_registration.csv     <- one row per retina
%   MaskComparison/README.txt
%
% Re-running the same key REPLACES that retina's row and images rather
% than duplicating them, so it is safe to redo one that looked wrong.
%
% Colours: RED his only, BLUE yours only, WHITE both.
%
% params:
%   workScale   fraction of his resolution to work at. Default 0.5, which
%               is ~1.38 um/pixel - still 2.7x finer than the QuRVA masks,
%               and a quarter of the memory. Use 1 for full resolution on
%               a retina you want publication-quality crops from.
%   pixSzeUM    um per pixel of HIS masks at full size. Default 1/1.4491
%   mode        'translation' | 'rigid' | 'similarity'. Default 'similarity'
%   workPx      long edge used for the registration estimate. Default 2000
%   overlayPx   long edge of the saved overlay. Default 2500
%   cropUM      crop width in microns. Default 300
%   exportMine  write your registered mask as TIFF. Default true

    if nargin < 4, params = struct(); end
    p.workScale  = getf(params,'workScale',  0.5);
    p.pixSzeUM   = getf(params,'pixSzeUM',   1/1.4491);
    p.mode       = getf(params,'mode',       'similarity');
    p.workPx     = getf(params,'workPx',     2000);
    p.overlayPx  = getf(params,'overlayPx',  2500);
    p.cropUM     = getf(params,'cropUM',     300);
    p.exportMine = getf(params,'exportMine', true);

    key = upper(strtrim(key));
    vascDir = fullfile(masterFolder,'VasculatureNumbers');
    maskDir = fullfile(masterFolder,'Masks');
    outRoot = fullfile(masterFolder,'MaskComparison');
    csvPath = fullfile(outRoot,'mask_registration.csv');

    for d = {'overlays','crops','my_masks'}
        dd = fullfile(outRoot, d{1});
        if ~exist(dd,'dir'), mkdir(dd); end
    end

    % ---------- find his file ----------
    hf = dir(fullfile(hisDir,'*whole_vessel_mask*.tif'));
    assert(~isempty(hf), 'No *whole_vessel_mask*.tif in %s', hisDir);
    hIdx = [];
    for k = 1:numel(hf)
        t = regexp(hf(k).name,'^([A-Za-z0-9]+)_(\d+)_whole_vessel','tokens','once');
        if ~isempty(t) && strcmp(upper([t{1} '_' t{2}]), key), hIdx = k; break; end
    end
    assert(~isempty(hIdx), 'No whole-vessel mask for %s in %s', key, hisDir);

    % ---------- find mine ----------
    mf = dir(fullfile(vascDir,'*.mat'));
    mIdx = [];
    for k = 1:numel(mf)
        t = regexp(mf(k).name,'^([A-Za-z0-9]+)_(\d+)','tokens','once');
        if ~isempty(t) && strcmp(upper([t{1} '_' t{2}]), key), mIdx = k; break; end
    end
    assert(~isempty(mIdx), 'No QuRVA .mat for %s in %s', key, vascDir);

    base = erase(mf(mIdx).name,'.mat');
    fprintf('\n%s\n  mine: %s\n  his : %s\n', key, base, hf(hIdx).name);

    % ---------- load ----------
    S = load(fullfile(vascDir, mf(mIdx).name),'smoothVessels');
    M = load(fullfile(maskDir, mf(mIdx).name),'thisMask');
    A0 = logical(S.smoothVessels); hull = logical(M.thisMask);
    clear S M

    B = imread(fullfile(hf(hIdx).folder, hf(hIdx).name));
    if size(B,3) > 1, B = any(B>0,3); else, B = B > 0; end

    if p.workScale < 1
        B = imresize(B, p.workScale, 'nearest');
    end
    pxUM = p.pixSzeUM / p.workScale;
    fprintf('  working at %d x %d, %.3f um/pixel\n', size(B,1), size(B,2), pxUM);

    A     = imresize(A0,   size(B),'nearest');
    hullB = imresize(hull, size(B),'nearest');
    clear A0 hull

    % ---------- register ----------
    s   = min(1, p.workPx / max(size(B)));
    tsm = imregcorr(imgaussfilt(double(imresize(A,s,'nearest')),2), ...
                    imgaussfilt(double(imresize(B,s,'nearest')),2), p.mode);
    Sd  = affine2d([s 0 0; 0 s 0; 0 0 1]);
    Su  = affine2d([1/s 0 0; 0 1/s 0; 0 0 1]);
    tform = affine2d(Sd.T * tsm.T * Su.T);

    ref  = imref2d(size(B));
    Areg = imwarp(A,     tform,'OutputView',ref,'InterpolationMethod','nearest');
    Hreg = imwarp(hullB, tform,'OutputView',ref,'InterpolationMethod','nearest');

    inside = Hreg & hullB;
    d0 = dicev(A & inside, B & inside);
    d1 = dicev(Areg & inside, B & inside);
    clear A hullB inside

    dx = tform.T(3,1); dy = tform.T(3,2);
    sc = hypot(tform.T(1,1), tform.T(1,2));
    rt = atan2d(tform.T(1,2), tform.T(1,1));

    Ain = Areg & Hreg;  Bin = B & Hreg;

    fprintf('  dice %.3f -> %.3f   shift %.1f um   scale %.4f   rot %+.2f deg\n', ...
            d0, d1, hypot(dx,dy)*pxUM, sc, rt);
    if d1 < 0.6
        warning('Dice is still %.2f after registration - check the overlay.', d1);
    end
    if abs(sc-1) > 0.02
        warning('Scale %.4f differs from 1 by more than 2%%.', sc);
    end

    % ---------- row ----------
    r = table(string(key), string(extractBefore(key,'_')), string(base), ...
              string(hf(hIdx).name), size(B,1), size(B,2), pxUM, p.workScale, ...
              d0, d1, d1-d0, dx, dy, dx*pxUM, dy*pxUM, hypot(dx,dy)*pxUM, ...
              sc, rt, string(p.mode), ...
              nnz(Ain)/max(nnz(Hreg),1), nnz(Bin)/max(nnz(Hreg),1), ...
              nnz(Ain&Bin), nnz(Bin&~Ain), nnz(Ain&~Bin), ...
        'VariableNames', {'key','strain','myFile','hisFile','imgHeight', ...
        'imgWidth','pxSzeUM','workScale','dice_before','dice_after', ...
        'dice_gain','shift_x_px','shift_y_px','shift_x_um','shift_y_um', ...
        'shift_mag_um','scale','rotation_deg','mode','myVesselFrac', ...
        'hisVesselFrac','bothPix','hisOnlyPix','myOnlyPix'});

    if exist(csvPath,'file')
        T = readtable(csvPath, 'TextType','string');
        T(strcmp(T.key, key), :) = [];      % replace, do not duplicate
        T = [T; r];
    else
        T = r;
    end
    T = sortrows(T, {'strain','key'});
    writetable(T, csvPath);

    % ---------- images ----------
    so = min(1, p.overlayPx / max(size(B)));
    imwrite(colr(imresize(Areg,so,'nearest'), imresize(Bin,so,'nearest')), ...
            fullfile(outRoot,'overlays', sprintf('%s_overlay.png', key)));

    z  = max(round(p.cropUM / pxUM), 64);
    st = regionprops(Hreg,'Centroid','BoundingBox');
    if ~isempty(st)
        ctr = st(1).Centroid; rad = 0.8*min(st(1).BoundingBox(3:4))/2;
        ang = [0 pi/2 pi 3*pi/2]; fr = [0.2 0.45 0.7 0.85];
        for k = 1:4
            cy = ctr(2)+rad*fr(k)*sin(ang(k));
            cx = ctr(1)+rad*fr(k)*cos(ang(k));
            r0 = max(min(round(cy-z/2), size(B,1)-z),1);
            c0 = max(min(round(cx-z/2), size(B,2)-z),1);
            imwrite(colr(Areg(r0:r0+z-1, c0:c0+z-1), Bin(r0:r0+z-1, c0:c0+z-1)), ...
                fullfile(outRoot,'crops', sprintf('%s_crop%d.png', key, k)));
        end
    end

    if p.exportMine
        imwrite(uint8(Areg)*255, fullfile(outRoot,'my_masks', ...
                sprintf('%s_qurva_mask_registered.tif', key)));
    end

    writeReadme(outRoot, p);

    fprintf('  %d retina(s) now in %s\n\n', height(T), csvPath);
end

% ======================================================================
function d = dicev(A,B), d = 2*nnz(A&B)/max(nnz(A)+nnz(B),1); end

function rgb = colr(A,B)
    rgb = cat(3, uint8(B)*255, uint8(A&B)*255, uint8(A)*255);
end

function writeReadme(outRoot, p)
    fid = fopen(fullfile(outRoot,'README.txt'),'w');
    fprintf(fid, [ ...
    'Mask comparison: QuRVA vs whole-retina masks\n' ...
    '============================================\n\n' ...
    'Last updated %s\n\n' ...
    'CONTENTS\n' ...
    '  overlays/    registered whole-retina overlay per retina\n' ...
    '  crops/       four crops per retina, centre to periphery\n' ...
    '  my_masks/    QuRVA masks registered onto the whole-retina grid\n' ...
    '  mask_registration.csv\n\n' ...
    'COLOURS\n' ...
    '  RED    whole-retina mask only\n' ...
    '  BLUE   QuRVA mask only\n' ...
    '  WHITE  both\n\n' ...
    'ALIGNMENT\n' ...
    '  The masks were not registered as received. Alignment was estimated\n' ...
    '  per retina by phase correlation (mode: %s) on the vessel pattern,\n' ...
    '  then applied to the QuRVA mask. dice_before and dice_after in the\n' ...
    '  CSV give the overlap either side of that correction.\n\n' ...
    'RESOLUTION\n' ...
    '  Comparison was done at %.0f%% of the whole-retina mask resolution\n' ...
    '  (%.3f um/pixel). QuRVA masks are coarser and were upsampled onto\n' ...
    '  that grid, so their edges are blocky - that is resampling, not\n' ...
    '  segmentation.\n\n' ...
    'READING THE OVERLAYS\n' ...
    '  Red forming complete vessels in black space = capillaries absent\n' ...
    '  from the QuRVA mask.\n' ...
    '  Red as a thin fringe along white vessels = boundary disagreement.\n' ...
    '  Two red lines through one white band = two capillaries merged into\n' ...
    '  one in the QuRVA mask.\n'], ...
    datestr(now), p.mode, 100*p.workScale, p.pixSzeUM/p.workScale);
    fclose(fid);
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end