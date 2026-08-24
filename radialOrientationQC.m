function radialOrientationQC(masterFolder, fileName, params)
% radialOrientationQC  Visual check of what radialOrientation is measuring.
%
%   radialOrientationQC(masterFolder, 'ACI_1_cage_12_Lectin(488).tif.mat')
%
% Produces two figures:
%
%  FIG 1 - branches coloured by their angle to the local outward direction
%          RED    = running radially (outward from the disc)
%          BLUE   = running circumferentially (around the disc)
%          GREY   = neither / 45 degrees
%          White circles are the radial bin boundaries, white * the apex.
%          If a visibly stringy outward region is not red, the measure is
%          not seeing what you see.
%
%  FIG 2 - zoom on one bin with the FITTED AXIS drawn on each branch as a
%          short white line, plus the local outward direction as a short
%          yellow line. These two lines being parallel means delta ~ 0
%          (radial); perpendicular means delta ~ 90 (circumferential).
%          This is the check that the axis estimation itself is right.

    if nargin < 3, params = struct(); end
    p.basePxUM       = getf(params,'basePxUM',       3/1.4491);
    p.nBins          = getf(params,'nBins',          6);
    p.radiusThreshUM = getf(params,'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params,'minBranchPix',   8);
    p.minLenUM       = getf(params,'minLenUM',       20);
    p.zoomBin        = getf(params,'zoomBin',        5);

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

    D = bwdist(~im_clean);
    skelRadius = zeros(size(D));
    skelRadius(im_skel) = D(im_skel);
    skel_cap = bwareaopen(im_skel & skelRadius > 0 & ...
                          skelRadius <= p.radiusThreshUM/pxUM, p.minBranchPix, 8);

    nbr    = conv2(double(skel_cap), [1 1 1;1 0 1;1 1 1], 'same');
    juncPx = skel_cap & nbr >= 3;
    L      = bwlabel(skel_cap & ~juncPx, 8);
    st     = regionprops(L, 'PixelList','Centroid','PixelIdxList');

    minPix = max(round(p.minLenUM / pxUM), 5);

    % ---- per-branch angle to the local outward direction ----
    [nr, nc] = size(im_hull);
    val = nan(numel(st),1);  cen = nan(numel(st),2);  ax = nan(numel(st),2);

    for k = 1:numel(st)
        P = st(k).PixelList;
        if size(P,1) < minPix, continue; end
        c = st(k).Centroid;
        Cv = cov(double(P));
        [V,Dg] = eig(Cv); [~,im_] = max(diag(Dg));
        v = V(:,im_); v = v/norm(v);
        r = [c(1)-apexXY(1); c(2)-apexXY(2)];
        if norm(r) < eps, continue; end
        r = r/norm(r);
        val(k) = 2*(abs(v'*r))^2 - 1;      % cos(2 delta): +1 radial, -1 circumf
        cen(k,:) = c;  ax(k,:) = v';
    end

    % ---- FIG 1: colour every branch by cos(2 delta) ----
    R = zeros(nr,nc); G = R; B = R;
    for k = 1:numel(st)
        if isnan(val(k)), continue; end
        t = (val(k)+1)/2;                   % 0 = circumferential, 1 = radial
        idx = st(k).PixelIdxList;
        R(idx) = 0.25 + 0.75*t;
        G(idx) = 0.35;
        B(idx) = 0.25 + 0.75*(1-t);
    end
    rgb = cat(3,R,G,B);
    rgb = max(rgb, 0.12*double(cat(3,im_clean,im_clean,im_clean)));

    [X,Y] = meshgrid(1:nc,1:nr);
    Dap   = hypot(X-apexXY(1), Y-apexXY(2));
    dHull = Dap(im_hull);
    edges = linspace(min(dHull), max(dHull)+eps, p.nBins+1);

    figure('Color','w','Name','radial angle per branch'); imshow(rgb); hold on
    th = linspace(0,2*pi,500);
    for b = 2:p.nBins
        plot(apexXY(1)+edges(b)*cos(th), apexXY(2)+edges(b)*sin(th), 'w-','LineWidth',1);
    end
    plot(apexXY(1), apexXY(2), 'w*','MarkerSize',12,'LineWidth',1.2);
    for b = 1:p.nBins
        band = im_hull & Dap>=edges(b) & Dap<edges(b+1);
        ro = radialOrientation(skel_cap, apexXY, band, pxUM, p.minLenUM);
        rm = (edges(b)+edges(b+1))/2;
        text(apexXY(1)+rm, apexXY(2), sprintf('%d\n%.2f', b, ro.radialAlignment), ...
             'Color','w','FontSize',11,'FontWeight','bold', ...
             'HorizontalAlignment','center');
    end
    title(sprintf('%s   red = radial, blue = circumferential', ...
          strrep(erase(fileName,'.mat'),'_','\_')));

    % ---- FIG 2: fitted axis vs outward direction, one bin ----
    b  = min(max(p.zoomBin,1), p.nBins);
    lo = edges(b); hi = edges(b+1);
    inB = ~isnan(val) & hypot(cen(:,1)-apexXY(1), cen(:,2)-apexXY(2)) >= lo & ...
          hypot(cen(:,1)-apexXY(1), cen(:,2)-apexXY(2)) < hi;
    ii = find(inB);
    if numel(ii) > 250, ii = ii(randperm(numel(ii),250)); end

    figure('Color','w','Name','fitted axis vs outward direction');
    imshow(0.18*double(cat(3,skel_cap,skel_cap,skel_cap))+0.05); hold on
    Lseg = 26;
    for k = ii(:)'
        c = cen(k,:); v = ax(k,:);
        plot(c(1)+[-1 1]*Lseg*v(1), c(2)+[-1 1]*Lseg*v(2), 'w-','LineWidth',1.4);
        r = [c(1)-apexXY(1), c(2)-apexXY(2)]; r = r/norm(r);
        plot(c(1)+[0 1]*Lseg*r(1), c(2)+[0 1]*Lseg*r(2), 'y-','LineWidth',1.1);
    end
    plot(apexXY(1)+edges(b)*cos(th),   apexXY(2)+edges(b)*sin(th),  'c-','LineWidth',1);
    plot(apexXY(1)+edges(b+1)*cos(th), apexXY(2)+edges(b+1)*sin(th),'c-','LineWidth',1);
    plot(apexXY(1), apexXY(2), 'w*','MarkerSize',12,'LineWidth',1.2);
    xlim(apexXY(1)+[-hi hi]*1.05); ylim(apexXY(2)+[-hi hi]*1.05);
    title(sprintf('bin %d: white = fitted branch axis, yellow = outward direction', b));

    fprintf('\n%s  px %.3f um  apex %.0f %.0f\n', erase(fileName,'.mat'), ...
            pxUM, apexXY(1), apexXY(2));
    for bb = 1:p.nBins
        band = im_hull & Dap>=edges(bb) & Dap<edges(bb+1);
        ro = radialOrientation(skel_cap, apexXY, band, pxUM, p.minLenUM);
        fprintf('  bin %d  %6.0f-%6.0f um   radAlign %+6.3f   globCoher %.3f   n %5d\n', ...
            bb, edges(bb)*pxUM, edges(bb+1)*pxUM, ...
            ro.radialAlignment, ro.globalCoherence, ro.nUsed);
    end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end