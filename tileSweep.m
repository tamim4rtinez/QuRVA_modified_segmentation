function T = tileSweep(masterFolder, params)
% tileSweep  How local coherence depends on window size, across strains.
%
%   T = tileSweep(newMaster, struct('basePxUM', 2/1.4491, 'radiusThreshUM', 9))
%
% Runs localCoherence at several tile sizes on a subset of retinas and
% writes tile_sweep.csv, one row per retina per tile size per bin.
%
% TWO THINGS TO READ OFF IT
%   medianBranches vs tileUM   - how many segments a window contains.
%       Below ~15 the coherence of a window is substantially chance-driven:
%       a handful of segments look parallel far more often than a crowd do.
%       This is the criterion to CHOOSE the tile size on, because it does
%       not involve the strain comparison.
%
%   meanCoh by bin, per tile   - how the radial pattern changes with window
%       size. Useful to see, but choosing the tile that maximises the strain
%       difference means selecting a parameter on the outcome. If you do
%       that, say so in the methods.
%
% params:
%   basePxUM        REQUIRED, um per pixel of the stored masks
%   tiles           tile widths in um. Default [150 200 250 300 400]
%   perStrain       retinas per strain. Default 2
%   radiusThreshUM  default 9
%   minBranchPix    default round(30/basePxUM)
%   minLenUM        default 20
%   minBranches     default 6
%   cohThresh       default 0.60

    if nargin < 2, params = struct(); end
    assert(isfield(params,'basePxUM'), 'params.basePxUM is required');
    p.basePxUM       = params.basePxUM;
    p.tiles          = getf(params,'tiles',          [150 200 250 300 400]);
    p.perStrain      = getf(params,'perStrain',      2);
    p.radiusThreshUM = getf(params,'radiusThreshUM', 9);
    p.minBranchPix   = getf(params,'minBranchPix',   round(30/params.basePxUM));
    p.minLenUM       = getf(params,'minLenUM',       20);
    p.minBranches    = getf(params,'minBranches',    6);
    p.cohThresh      = getf(params,'cohThresh',      0.60);

    f = dir(fullfile(masterFolder,'VasculatureNumbers','*.mat'));
    names   = {f.name};
    strains = cellfun(@(s) regexprep(s,'_.*$',''), names, 'UniformOutput', false);
    uS      = unique(strains,'stable');

    picks = {};
    for k = 1:numel(uS)
        idx = find(strcmp(strains, uS{k}));
        picks = [picks, names(idx(1:min(p.perStrain, numel(idx))))]; %#ok<AGROW>
    end
    fprintf('\n%d retinas x %d tile sizes = %d runs\n\n', ...
            numel(picks), numel(p.tiles), numel(picks)*numel(p.tiles));

    rows = [];
    for k = 1:numel(picks)
        for t = p.tiles
            lp = struct('basePxUM', p.basePxUM, 'radiusThreshUM', p.radiusThreshUM, ...
                        'minBranchPix', p.minBranchPix, 'minLenUM', p.minLenUM, ...
                        'tileUM', t, 'minBranches', p.minBranches, ...
                        'cohThresh', p.cohThresh, 'showQC', false);
            try
                [Tb, m] = localCoherence(masterFolder, picks{k}, lp);
                for b = 1:height(Tb)
                    r.fileName       = string(erase(picks{k},'.mat'));
                    r.strain         = string(regexprep(picks{k},'_.*$',''));
                    r.tileUM         = t;
                    r.bin            = Tb.bin(b);
                    r.nTiles         = Tb.nTiles(b);
                    r.meanCoh        = Tb.meanCoh(b);
                    r.p90Coh         = Tb.p90Coh(b);
                    r.fracAligned    = Tb.fracAligned(b);
                    r.medianBranches = median(double(m.nbr));
                    r.totalWindows   = numel(m.coh);
                    rows = [rows; r]; %#ok<AGROW>
                end
                fprintf('  %-22s tile %3d um: %4d windows, median %2d branches\n', ...
                        erase(picks{k},'.mat'), t, numel(m.coh), median(m.nbr));
            catch ME
                fprintf('  %-22s tile %3d um FAILED: %s\n', picks{k}, t, ME.message);
            end
        end
    end

    T = struct2table(rows);
    writetable(T, fullfile(masterFolder,'tile_sweep.csv'));
    fprintf('\nWrote %s\n', fullfile(masterFolder,'tile_sweep.csv'));

    % ---- branches per window vs tile size, one line per strain ----
    G = groupsummary(T(T.bin==1,:), {'strain','tileUM'}, 'mean', 'medianBranches');
    figure('Color','w','Position',[80 80 1100 420]);

    subplot(1,2,1); hold on
    uu = unique(G.strain);
    for k = 1:numel(uu)
        g = G(G.strain==uu(k),:);
        plot(g.tileUM, g.mean_medianBranches, '-o','LineWidth',1.2);
    end
    yline(15,'--','15 branches','Color',[.4 .4 .4]);
    xlabel('tile width (\mum)'); ylabel('median branches per window');
    legend(uu,'Location','northwest'); box off
    title('Window occupancy');

    % ---- peripheral minus central contrast vs tile size ----
    subplot(1,2,2); hold on
    for k = 1:numel(uu)
        c = zeros(numel(p.tiles),1);
        for j = 1:numel(p.tiles)
            s = T(T.strain==uu(k) & T.tileUM==p.tiles(j), :);
            c(j) = mean(s.meanCoh(s.bin>=5)) - mean(s.meanCoh(s.bin<=3));
        end
        plot(p.tiles, c, '-o','LineWidth',1.2);
    end
    yline(0,'-','Color',[.7 .7 .7]);
    xlabel('tile width (\mum)'); ylabel('meanCoh: bins 5-6 minus bins 1-3');
    legend(uu,'Location','northwest'); box off
    title('Radial contrast');
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end