function maskOverlayQC(masterFolder, otherMaskPath, fileName, params)
% maskOverlayQC  Colour-coded overlay of two vessel masks. Look before
%                you measure.
%
%   maskOverlayQC(masterFolder, '/path/to/his/ACI_1.tif', ...
%                 'ACI_1_cage_12_Lectin(488).tif.mat')
%
% COLOURS
%   RED     his mask only
%   BLUE    your (QuRVA) mask only
%   WHITE   both
%   black   neither
%
% Writes to masterFolder/MaskComparison/<base>/ :
%   whole_retina.png   downscaled, for orientation only
%   crop_1..N.png      FULL RESOLUTION crops - these are the ones to judge from
%   crop_map.png       whole retina with the crop boxes marked
%
% WHAT TO LOOK FOR
%   Red forming complete vessels running through black space
%       -> genuine capillaries your mask missed.
%   Red as a thin fringe along white vessels
%       -> boundary disagreement only, not missed vessels.
%   Two red lines running parallel through one white band
%       -> your mask merged two capillaries into one. This is the failure
%          mode that would shorten your skeleton without losing area.
%   Red and blue running parallel along every vessel
%       -> the masks are not registered. Stop and fix alignment first.
%
% params:
%   nCrops      number of full-resolution crops. Default 6
%   cropUM      crop width in microns. Default 300
%   pixSzeUM    um per pixel of HIS mask. Default 1/1.4491 (0.690)
%   wholeMaxPx  long edge of the downscaled overview. Default 2000
%   cropAt      Nx2 [x y] crop centres in HIS pixel coords. Default: spread
%               across the retina at increasing eccentricity

    if nargin < 4, params = struct(); end
    p.nCrops     = getf(params,'nCrops',     6);
    p.cropUM     = getf(params,'cropUM',     300);
    p.pixSzeUM   = getf(params,'pixSzeUM',   1/1.4491);
    p.wholeMaxPx = getf(params,'wholeMaxPx', 2000);
    p.cropAt     = getf(params,'cropAt',     []);
        p.tform = getf(params,'tform', []);

    S = load(fullfile(masterFolder,'VasculatureNumbers',fileName), 'smoothVessels');
    M = load(fullfile(masterFolder,'Masks',fileName), 'thisMask');
    A0   = logical(S.smoothVessels);      % yours, QuRVA resolution
    hull = logical(M.thisMask);

    B = imread(otherMaskPath);            % his, full resolution
    if size(B,3) > 1, B = any(B > 0, 3); else, B = B > 0; end

    fprintf('\nyours %d x %d   his %d x %d   ratio %.2f\n', ...
        size(A0,1), size(A0,2), size(B,1), size(B,2), size(B,1)/size(A0,1));

    % put yours on his grid - never the other way round, downsampling his
    % would erase exactly the thin vessels we are looking for
      A    = imresize(A0,   size(B), 'nearest');
    hull = imresize(hull, size(B), 'nearest');

    if ~isempty(p.tform)
        ref  = imref2d(size(B));
        A    = imwarp(A,    p.tform, 'OutputView', ref, 'InterpolationMethod','nearest');
        hull = imwarp(hull, p.tform, 'OutputView', ref, 'InterpolationMethod','nearest');
        fprintf('Registration applied.\n');
    end

    A    = A & hull;
    B    = B & hull;

    dice = 2*nnz(A & B) / max(nnz(A)+nnz(B),1);
       fprintf('Dice %.3f%s   yours %.2f%% of retina   his %.2f%% of retina\n', ...
        dice, ternary(isempty(p.tform),'',' (registered)'), ...
        100*nnz(A)/nnz(hull), 100*nnz(B)/nnz(hull));
    if dice < 0.6
        warning(['Dice %.2f is low - the masks may not be aligned. ' ...
                 'Check the whole-retina image before trusting the crops.'], dice);
    end

    base   = erase(fileName,'.mat');
    outDir = fullfile(masterFolder,'MaskComparison', matlab.lang.makeValidName(base));
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % ---- whole retina, downscaled ----
    s = min(1, p.wholeMaxPx / max(size(B)));
    rgbSmall = colourise(imresize(A,s,'nearest'), imresize(B,s,'nearest'));
    imwrite(rgbSmall, fullfile(outDir,'whole_retina.png'));

    % ---- choose crop centres ----
    z = max(round(p.cropUM / p.pixSzeUM), 64);
    if isempty(p.cropAt)
        st  = regionprops(hull,'Centroid','BoundingBox');
        ctr = st(1).Centroid;
        rad = 0.85 * min(st(1).BoundingBox(3:4)) / 2;
        ang = linspace(0, 2*pi, p.nCrops+1); ang(end) = [];
        frac = linspace(0.15, 0.85, p.nCrops);
        cx = ctr(1) + rad*frac(:).*cos(ang(:));
        cy = ctr(2) + rad*frac(:).*sin(ang(:));
        p.cropAt = [cx cy];
    end

    boxes = zeros(size(p.cropAt,1), 4);
    for k = 1:size(p.cropAt,1)
        r0 = max(min(round(p.cropAt(k,2)-z/2), size(B,1)-z), 1);
        c0 = max(min(round(p.cropAt(k,1)-z/2), size(B,2)-z), 1);
        rr = r0:r0+z-1; cc = c0:c0+z-1;
        boxes(k,:) = [c0 r0 z z];

        rgb = colourise(A(rr,cc), B(rr,cc));
        imwrite(rgb, fullfile(outDir, sprintf('crop_%d.png', k)));

        both = nnz(A(rr,cc) & B(rr,cc));
        hOnly = nnz(B(rr,cc) & ~A(rr,cc));
        yOnly = nnz(A(rr,cc) & ~B(rr,cc));
        fprintf('  crop %d at [%5.0f %5.0f]  both %6d  his-only %6d  yours-only %6d\n', ...
                k, p.cropAt(k,1), p.cropAt(k,2), both, hOnly, yOnly);
    end

    % ---- map showing where the crops are ----
    f = figure('Visible','off','Color','w');
    imshow(rgbSmall); hold on
    for k = 1:size(boxes,1)
        rectangle('Position', boxes(k,:)*s, 'EdgeColor','y', 'LineWidth',1);
        text(boxes(k,1)*s, boxes(k,2)*s - 8, sprintf('%d',k), ...
             'Color','y','FontWeight','bold');
    end
    exportgraphics(f, fullfile(outDir,'crop_map.png'), 'Resolution',150);
    close(f);

    fprintf(['\nRED = his only, BLUE = yours only, WHITE = both.\n' ...
             'Wrote %s\n\n'], outDir);
end

% ======================================================================
function rgb = colourise(A, B)
    R = uint8(B)*255;                 % his
    Bl= uint8(A)*255;                 % yours
    G = uint8(A & B)*255;             % overlap -> adds green, giving white
    rgb = cat(3, R, G, Bl);
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function out = ternary(c,a,b), if c, out = a; else, out = b; end, end