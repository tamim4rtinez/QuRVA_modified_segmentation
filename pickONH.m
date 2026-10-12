function pickONH(masterFolder, params)
% pickONH  Click the optic nerve head on each retina and save the position.
%
%   pickONH(newMaster)
%   pickONH(newMaster, struct('redo', true))
%
% Writes masterFolder/ONCenter/<name>.mat containing thisONCenter = [x y],
% in the same pixel coordinates as the mask - which is what every radial
% script expects.
%
% CONTROLS
%   click        set the optic nerve head, move to the next retina
%   right-click  skip this retina for now
%   z            undo, go back one retina
%   q            stop (progress is saved as you go, so you can resume)
%
% Retinas that already have a saved position are skipped unless redo is set.
%
% params:
%   redo        revisit retinas that already have a position. Default false
%   dispPx      long edge of the displayed image. Default 1400. The saved
%               coordinate is converted back to full mask resolution.
%   pattern     which retinas to do, e.g. 'ACI*'. Default '*'

    if nargin < 2, params = struct(); end
    p.redo    = getf(params,'redo',    false);
    p.dispPx  = getf(params,'dispPx',  1400);
    p.pattern = getf(params,'pattern', '*');

    vascDir = fullfile(masterFolder,'VasculatureNumbers');
    ocDir   = fullfile(masterFolder,'ONCenter');
    assert(exist(vascDir,'dir')==7, 'No VasculatureNumbers in %s', masterFolder);
    if ~exist(ocDir,'dir'), mkdir(ocDir); end

    files = dir(fullfile(vascDir, [p.pattern '.mat']));
    assert(~isempty(files), 'No .mat files matching %s', p.pattern);

    todo = true(numel(files),1);
    if ~p.redo
        for k = 1:numel(files)
            todo(k) = ~exist(fullfile(ocDir, files(k).name), 'file');
        end
    end
    idx = find(todo);
    fprintf('\n%d retina(s) to do (%d already have a position)\n', ...
            numel(idx), numel(files)-numel(idx));
    if isempty(idx), return; end

    fprintf(['click = set and advance   right-click = skip   ' ...
             'z = back   q = quit\n\n']);

    fig = figure('Name','Click the optic nerve head','Color','w', ...
                 'NumberTitle','off');
    i = 1;
    while i >= 1 && i <= numel(idx)
        k = idx(i);
        name = files(k).name;

        S = load(fullfile(vascDir,name), 'smoothVessels');
        BW = S.smoothVessels;
        s  = min(1, p.dispPx / max(size(BW)));
        Bs = imresize(BW, s, 'nearest');

        clf(fig);
        imshow(Bs); hold on
        if exist(fullfile(ocDir,name),'file')
            C = load(fullfile(ocDir,name),'thisONCenter');
            plot(C.thisONCenter(1)*s, C.thisONCenter(2)*s, ...
                 'r+','MarkerSize',18,'LineWidth',2);
        end
        title(sprintf('%d / %d   %s', i, numel(idx), strrep(name,'_','\_')), ...
              'Interpreter','tex');

        [x, y, btn] = ginput(1);

        if isempty(btn)                      % window closed
            fprintf('Stopped at %d / %d.\n', i, numel(idx)); break
        elseif btn == 'q' || btn == 'Q'
            fprintf('Stopped at %d / %d.\n', i, numel(idx)); break
        elseif btn == 'z' || btn == 'Z'
            i = max(i-1, 1); continue
        elseif btn == 3                      % right-click, skip
            fprintf('  skipped  %s\n', name); i = i + 1; continue
        end

        thisONCenter = [x y] / s;            %#ok<NASGU>  back to mask pixels
        save(fullfile(ocDir, name), 'thisONCenter');
        fprintf('  %-30s  [%6.0f %6.0f]\n', erase(name,'.mat'), x/s, y/s);
        i = i + 1;
    end

    if ishandle(fig), close(fig); end

    done = numel(dir(fullfile(ocDir,'*.mat')));
    fprintf('\n%d of %d retinas now have an optic nerve head position.\n\n', ...
            done, numel(files));
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end