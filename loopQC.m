function loopQC(masterFolder, params)
% loopQC  Show exactly which enclosed regions are being counted as loops.
%
%   loopQC(masterFolder)
%   loopQC(masterFolder, struct('zoomUM', 900, 'minLoopUM2', 100))
%
% Takes the first retina of each strain and draws, for each:
%   grey    capillary skeleton
%   COLOUR  loops that PASS the size filter, tinted by area
%           (blue = small, yellow = large)
%   red     enclosed regions REJECTED as too small - these are the
%           diagonal-staircase artifacts. If red dominates the image,
%           the filter is doing the heavy lifting and the threshold
%           deserves scrutiny.
%
% A second figure zooms one square region per strain at full pixel
% resolution, so you can see whether a counted loop is a real mesh
% space or two skeleton branches touching.

    if nargin < 2, params = struct(); end
    p.basePxUM       = getf(params,'basePxUM',       3/1.4491);
    p.radiusThreshUM = getf(params,'radiusThreshUM', 17.1);
    p.minBranchPix   = getf(params,'minBranchPix',   8);
    p.minLoopUM2     = getf(params,'minLoopUM2',     100);
    p.zoomUM         = getf(params,'zoomUM',         900);
    p.gapPix = getf(params,'gapPix', 3);
    p.saveDir = getf(params,'saveDir', fullfile(masterFolder,'CapillaryQC'));
    if ~exist(p.saveDir,'dir'), mkdir(p.saveDir); end

    files   = dir(fullfile(masterFolder,'VasculatureNumbers','*.mat'));
    names   = {files.name};
    strains = cellfun(@(s) extractBefore(s,'_'), names, 'UniformOutput', false);
    uS      = unique(strains, 'stable');

   
    fprintf('\n%-8s %-40s %7s %7s %7s %9s %9s\n', 'strain','file', ...
            'kept','reject','%%kept','medArea','maxArea');
    fprintf('%s\n', repmat('-',1,90));

    for k = 1:numel(uS)
        f = names{find(strcmp(strains, uS{k}), 1)};

        S = load(fullfile(masterFolder,'VasculatureNumbers',f), ...
                 'smoothVessels','vesselSkelMask');
        M = load(fullfile(masterFolder,'Masks',f), 'thisMask');

        im_clean = logical(S.smoothVessels);
        im_skel  = logical(S.vesselSkelMask);
        im_hull  = logical(M.thisMask);

        sf = 1;
        if ~isequal(size(im_hull), size(im_clean))
            sf = size(im_hull,1)/size(im_clean,1);
            im_hull = imresize(im_hull, size(im_clean), 'nearest');
        end
        pxUM = p.basePxUM * sf;

        D  = bwdist(~im_clean);
        sr = zeros(size(D)); sr(im_skel) = D(im_skel);
        skel_cap = bwareaopen(im_skel & sr > 0 & ...
                              sr <= p.radiusThreshUM/pxUM, p.minBranchPix, 8);

        % ---- enclosed regions of the skeleton ----
        skel_fix = bwmorph(imclose(skel_cap, strel('disk', p.gapPix)), 'thin', Inf);
        allLoops = imfill(skel_fix,'holes') & ~skel_fix;
        minPix   = max(round(p.minLoopUM2 / pxUM^2), 2);

        cc  = bwconncomp(allLoops, 8);
        A   = cellfun(@numel, cc.PixelIdxList)';
        keep = A >= minPix;

        keptMask = false(size(allLoops));
        keptMask(vertcat(cc.PixelIdxList{keep}))  = true;
        rejMask  = false(size(allLoops));
        if any(~keep)
            rejMask(vertcat(cc.PixelIdxList{~keep})) = true;
        end

        areas_um2 = A(keep) * pxUM^2;

        fprintf('%-8s %-40s %7d %7d %6.1f%% %9.0f %9.0f\n', uS{k}, ...
            erase(f,'.mat'), nnz(keep), nnz(~keep), ...
            100*nnz(keep)/max(numel(A),1), ...
            median(areas_um2), max(areas_um2));

        % ---- tint kept loops by area ----
        Lk = zeros(size(allLoops));
        idxK = find(keep);
        for j = 1:numel(idxK)
            Lk(cc.PixelIdxList{idxK(j)}) = A(idxK(j));
        end
        tint = zeros(size(Lk));
        if any(Lk(:) > 0)
            tint(Lk>0) = min(log10(Lk(Lk>0)) / log10(prctile(A(keep),98)), 1);
        end

        R = 0.30*double(skel_fix); G = R; B = R;
        R(keptMask) = 0.15 + 0.85*tint(keptMask);
        G(keptMask) = 0.35 + 0.60*tint(keptMask);
        B(keptMask) = 0.95 - 0.75*tint(keptMask);
        R(rejMask)  = 0.95; G(rejMask) = 0.10; B(rejMask) = 0.10;
        rgb = cat(3,R,G,B);
        fh = figure('Color','w','Visible','on','Position',[80 80 1400 1200]);
        imshow(rgb);
        title(sprintf('%s   %d loops kept, %d rejected (%.0f%% kept)', ...
              uS{k}, nnz(keep), nnz(~keep), ...
              100*nnz(keep)/max(numel(A),1)), 'FontSize', 14);
        exportgraphics(fh, fullfile(p.saveDir, ...
            sprintf('loopQC_%s.jpg', uS{k})), 'Resolution', 200);
        close(fh);
       
    end

    fprintf(['\nRed = rejected as smaller than %g um2. Colour = kept,\n' ...
             'blue small to yellow large. Check the zoom: a counted loop\n' ...
             'should enclose visible open space, not just two branches\n' ...
             'touching at a corner.\n\n'], p.minLoopUM2);
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end