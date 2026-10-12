function showRadialBins(masterFolder, fileName, nBins, highlight, params)
% showRadialBins  Draw the radial bin boundaries on one retina.
%
%   showRadialBins(masterFolder, 'ACI_3_cage2_Lectin(488).tif.mat')
%   showRadialBins(masterFolder, fname, 6, [3 4])
%
% Reproduces the binning used by capillaryMetrics/radialHeterogeneity so
% you can see exactly where each annulus falls. Capillary skeleton is
% green, large-vessel skeleton red, bin boundaries white, highlighted
% bins shaded. Prints the per-bin density so it matches the CSV.

    if nargin < 3 || isempty(nBins),     nBins = 6;        end
    if nargin < 4,                        highlight = [];   end

    S = load(fullfile(masterFolder,'VasculatureNumbers',fileName), ...
             'smoothVessels','vesselSkelMask');
    M = load(fullfile(masterFolder,'Masks',fileName), 'thisMask');
    C = load(fullfile(masterFolder,'ONCenter',fileName), 'thisONCenter');

    im_clean = logical(S.smoothVessels);
    im_skel  = logical(S.vesselSkelMask);
    im_hull  = logical(M.thisMask);

    % Match mask resolution to vessel mask, exactly as statsFromQuRVA does
    scaleFactor = 1;
    if ~isequal(size(im_hull), size(im_clean))
        scaleFactor = size(im_hull,1) / size(im_clean,1);
        im_hull = imresize(im_hull, size(im_clean), 'nearest');
    end
    apexXY = C.thisONCenter(:).' / scaleFactor;

        if nargin < 5, params = struct(); end
    p.basePxUM       = getf(params,'basePxUM',       []);
    p.radiusThreshUM = getf(params,'radiusThreshUM', []);
    p.radiusThreshPix= getf(params,'radiusThreshPix',4.5);
    p.minBranchPix   = getf(params,'minBranchPix',   8);

    % physical threshold if given, else pixels
    if ~isempty(p.radiusThreshUM) && ~isempty(p.basePxUM)
        pxUM            = p.basePxUM * scaleFactor;
        radiusThreshPix = p.radiusThreshUM / pxUM;
    else
        pxUM            = NaN;
        radiusThreshPix = p.radiusThreshPix;
    end
    minBranchPix = p.minBranchPix;

    % Capillary / large split (same rule as capillaryMetrics)
    radiusThreshPix = 4.5;
    minBranchPix    = 8;
    D = bwdist(~im_clean);
    skelRadius = zeros(size(D));
    skelRadius(im_skel) = D(im_skel);
    skel_cap   = im_skel & skelRadius > 0 & skelRadius <= radiusThreshPix;
    skel_large = im_skel & skelRadius > radiusThreshPix;
    skel_cap   = bwareaopen(skel_cap, minBranchPix, 8);

    % Radial bins
    [rr, cc] = size(im_hull);
    [X, Y]   = meshgrid(1:cc, 1:rr);
    Dap      = hypot(X - apexXY(1), Y - apexXY(2));
    dHull    = Dap(im_hull);
    edges    = linspace(min(dHull), max(dHull)+eps, nBins+1);

    % ---- figure ----
    rgb = repmat(0.18*double(im_clean), [1 1 3]);
    rgb(:,:,2) = max(rgb(:,:,2), double(imdilate(skel_cap,   strel('disk',1))));
    rgb(:,:,1) = max(rgb(:,:,1), double(imdilate(skel_large, strel('disk',1))));
    rgb(:,:,3) = max(rgb(:,:,3), double(imdilate(skel_large, strel('disk',1))));

    figure('Color','w'); imshow(rgb); hold on
    
    function v = getf(s,f,d)
        if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
    end
    % shade highlighted bins
    for b = highlight(:)'
        band = im_hull & Dap >= edges(b) & Dap < edges(b+1);
        h = imshow(cat(3, ones(rr,cc), ones(rr,cc), zeros(rr,cc)));
        set(h, 'AlphaData', 0.16*double(band));
    end

    % outline the retina and draw the ring boundaries
    visboundaries(im_hull, 'Color','w', 'LineWidth',0.5, 'EnhanceVisibility',false);
    th = linspace(0, 2*pi, 400);
    for b = 2:nBins +1
        plot(apexXY(1)+edges(b)*cos(th), apexXY(2)+edges(b)*sin(th), ...
             'w-', 'LineWidth', 1);
    end
    plot(apexXY(1), apexXY(2), 'w*', 'MarkerSize', 10, 'LineWidth', 1.2);

    % bin numbers along the horizontal from the apex
    for b = 1:nBins
        rmid = (edges(b) + edges(b+1))/2;
        text(apexXY(1)+rmid, apexXY(2), sprintf('%d', b), ...
             'Color','w','FontSize',12,'FontWeight','bold', ...
             'HorizontalAlignment','center');
    end

    title(strrep(fileName,'_','\_'), 'Interpreter','tex');

    % ---- per-bin density, same numbers as the CSV ----
    fprintf('\n%-28s  %10s %10s %12s\n', fileName, 'skelPix', 'tissuePix', 'density');
    for b = 1:nBins
        band = im_hull & Dap >= edges(b) & Dap < edges(b+1);
        hp   = nnz(band);
        sp   = nnz(skel_cap & band);
        fprintf('  bin %d %22s %10d %10d %12.5f\n', b, '', sp, hp, sp/max(hp,1));
    end
end