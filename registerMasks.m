function [tform, info] = registerMasks(masterFolder, otherMaskPath, fileName, params)
% registerMasks  Find the transform that aligns QuRVA's mask to another
%                pipeline's mask, and report how much it helps.
%
%   [tform, info] = registerMasks(masterFolder, hisMaskPath, fileName)
%
% Estimates the alignment by phase correlation rather than by landmarks.
% A single landmark (e.g. the optic nerve head) can only correct
% translation - if the two masks also differ slightly in scale or rotation,
% aligning one point leaves the periphery misaligned. Phase correlation
% uses the whole vessel pattern, so it catches all three.
%
% Runs on downscaled copies for speed, then rescales the transform to full
% resolution. A 0.69 um/pixel whole-retina mask is ~434 Mpx; registering at
% that size is unnecessary and slow.
%
% Writes before/after overlay crops to
%   masterFolder/MaskComparison/<base>/registration/
%
% params:
%   mode        'translation' | 'rigid' | 'similarity'. Default 'similarity'
%               (translation + rotation + uniform scale). Use 'translation'
%               if you are confident the scale is already right.
%   workPx      long edge used for estimation. Default 2000
%   cropUM      width of the check crops. Default 300
%   pixSzeUM    um per pixel of HIS mask. Default 1/1.4491
%
% RETURNS
%   tform   affine2d mapping YOUR mask onto HIS, at his resolution
%   info    struct with dice before and after, and the shift in um
%
% Feed tform into the comparison step with imwarp:
%   Areg = imwarp(A, tform, 'OutputView', imref2d(size(B)), ...
%                 'InterpolationMethod','nearest');

    if nargin < 4, params = struct(); end
    p.mode     = getf(params,'mode',     'similarity');
    p.workPx   = getf(params,'workPx',   2000);
    p.cropUM   = getf(params,'cropUM',   300);
    p.pixSzeUM = getf(params,'pixSzeUM', 1/1.4491);

    S = load(fullfile(masterFolder,'VasculatureNumbers',fileName),'smoothVessels');
    M = load(fullfile(masterFolder,'Masks',fileName),'thisMask');
    A0   = logical(S.smoothVessels);
    hull = logical(M.thisMask);

    B = imread(otherMaskPath);
    if size(B,3) > 1, B = any(B > 0, 3); else, B = B > 0; end

    A = imresize(A0, size(B), 'nearest');      % yours onto his grid
    H = imresize(hull, size(B), 'nearest');

    fprintf('\n%s\n', fileName);
    fprintf('  grid %d x %d\n', size(B,1), size(B,2));

    % ---- estimate on downscaled copies ----
    s  = min(1, p.workPx / max(size(B)));
    As = imresize(A, s, 'nearest');
    Bs = imresize(B, s, 'nearest');

    % light blur helps phase correlation lock onto structure rather than
    % individual pixel edges
    Af = imgaussfilt(double(As), 2);
    Bf = imgaussfilt(double(Bs), 2);

    tsmall = imregcorr(Af, Bf, p.mode);

    % rescale the transform from the working grid back to full resolution
    Sdown = affine2d([s 0 0; 0 s 0; 0 0 1]);
    Sup   = affine2d([1/s 0 0; 0 1/s 0; 0 0 1]);
    tform = affine2d(Sdown.T * tsmall.T * Sup.T);

    % ---- how much did it help ----
    ref   = imref2d(size(B));
    Areg  = imwarp(A, tform, 'OutputView', ref, 'InterpolationMethod','nearest');
    Hreg  = imwarp(H, tform, 'OutputView', ref, 'InterpolationMethod','nearest');
    inside = Hreg & H;

    d0 = dice(A    & inside, B & inside);
    d1 = dice(Areg & inside, B & inside);

    dx = tform.T(3,1); dy = tform.T(3,2);
    sc = hypot(tform.T(1,1), tform.T(1,2));
    rot = atan2d(tform.T(1,2), tform.T(1,1));

    fprintf('  mode %s\n', p.mode);
    fprintf('  shift  %+.1f, %+.1f px  (%+.1f, %+.1f um)\n', ...
            dx, dy, dx*p.pixSzeUM, dy*p.pixSzeUM);
    fprintf('  scale  %.4f    rotation %+.2f deg\n', sc, rot);
    fprintf('  dice   %.3f -> %.3f\n', d0, d1);

    if d1 < d0 + 0.01
        warning(['Registration barely improved the overlap. Either the ' ...
                 'masks were already aligned, or the offset is not a ' ...
                 'global transform (e.g. a crop difference). Look at the ' ...
                 'before/after crops.']);
    end
    if abs(sc-1) > 0.02
        warning(['Estimated scale is %.4f, not 1. The two masks are at ' ...
                 'genuinely different scales - worth resolving before ' ...
                 'comparing anything.'], sc);
    end

    info = struct('dice_before',d0,'dice_after',d1,'dx_px',dx,'dy_px',dy, ...
                  'scale',sc,'rot_deg',rot,'mode',p.mode);

    % ---- before / after crops ----
    base   = erase(fileName,'.mat');
    outDir = fullfile(masterFolder,'MaskComparison', ...
                      matlab.lang.makeValidName(base),'registration');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    z  = max(round(p.cropUM / p.pixSzeUM), 64);
    st = regionprops(H,'Centroid','BoundingBox');
    ctr = st(1).Centroid;  rad = 0.8*min(st(1).BoundingBox(3:4))/2;
    ang = linspace(0,2*pi,5); ang(end) = [];
    fr  = [0.2 0.45 0.7 0.85];

    for k = 1:4
        cx = ctr(1) + rad*fr(k)*cos(ang(k));
        cy = ctr(2) + rad*fr(k)*sin(ang(k));
        r0 = max(min(round(cy-z/2), size(B,1)-z),1);
        c0 = max(min(round(cx-z/2), size(B,2)-z),1);
        rr = r0:r0+z-1; cc = c0:c0+z-1;
        imwrite(col(A(rr,cc),    B(rr,cc)), fullfile(outDir, sprintf('crop%d_before.png',k)));
        imwrite(col(Areg(rr,cc), B(rr,cc)), fullfile(outDir, sprintf('crop%d_after.png',k)));
    end

    fprintf('  RED = his only, BLUE = yours only, WHITE = both\n');
    fprintf('  Wrote %s\n\n', outDir);
end

% ======================================================================
function d = dice(A,B), d = 2*nnz(A&B)/max(nnz(A)+nnz(B),1); end

function rgb = col(A,B)
    rgb = cat(3, uint8(B)*255, uint8(A&B)*255, uint8(A)*255);
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end