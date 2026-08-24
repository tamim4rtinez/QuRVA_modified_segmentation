function [Tbins, Twin] = batchLocalCoherence(masterFolder, params)
% batchLocalCoherence  Run localCoherence over every retina.
%
%   [Tbins, Twin] = batchLocalCoherence(masterFolder)
%
% Writes two CSVs into masterFolder:
%
%   local_coherence_bins.csv     one row per retina per bin
%       meanCoh, p90Coh, fracAligned, nTiles, medianBranchesPerTile
%       -> use for the ANOVA and the bar / radial-profile plots
%
%   local_coherence_windows.csv  one row per window
%       fileName, strain, bin, coh, nBranches
%       -> use for distribution plots, and to check whether coherence
%          is confounded with how many branches a window contains
%
% NOTE: windows overlap by 50%, so window rows are NOT independent.
% Statistics go on the per-retina values in Tbins, never on window rows.

    if nargin < 2, params = struct(); end
    params.showQC = false;

    files = dir(fullfile(masterFolder,'VasculatureNumbers','*.mat'));
    assert(~isempty(files), 'No .mat files found');

    Tbins = [];  Twin = [];

    for ii = 1:numel(files)
        f = files(ii).name;
        try
            [T, m] = localCoherence(masterFolder, f, params);

            base   = erase(f, '.mat');
            strain = extractBefore(base, '_');
            T.strain = repmat(string(strain), height(T), 1);
            Tbins = [Tbins; T]; %#ok<AGROW>

            d = hypot(m.tx - m.apex(1), m.ty - m.apex(2));
            wbin = zeros(numel(d),1);
            for b = 1:numel(m.edges)-1
                wbin(d >= m.edges(b) & d < m.edges(b+1)) = b;
            end
            keep = wbin > 0;

            W = table(repmat(string(base),   nnz(keep), 1), ...
                      repmat(string(strain), nnz(keep), 1), ...
                      wbin(keep), m.coh(keep), m.nbr(keep), ...
                      'VariableNames', {'fileName','strain','bin','coh','nBranches'});
            Twin = [Twin; W]; %#ok<AGROW>

            fprintf('[%2d/%2d] %-42s  %d windows\n', ii, numel(files), base, nnz(keep));

        catch ME
            fprintf('[%2d/%2d] %-42s  FAILED: %s\n', ii, numel(files), f, ME.message);
        end
    end

    % median branches per tile, per retina per bin - for the confound check
    G = groupsummary(Twin, {'fileName','bin'}, 'median', 'nBranches');
    G.Properties.VariableNames{end} = 'medianBranchesPerTile';
    Tbins = outerjoin(Tbins, G(:, {'fileName','bin','medianBranchesPerTile'}), ...
                      'Keys', {'fileName','bin'}, 'MergeKeys', true, 'Type', 'left');

    writetable(Tbins, fullfile(masterFolder,'local_coherence_bins.csv'));
    writetable(Twin,  fullfile(masterFolder,'local_coherence_windows.csv'));

    fprintf('\nWrote local_coherence_bins.csv (%d rows)\n', height(Tbins));
    fprintf('Wrote local_coherence_windows.csv (%d rows)\n', height(Twin));

    fprintf('\n--- fracAligned by strain and bin ---\n');
    S = groupsummary(Tbins, {'strain','bin'}, 'mean', 'fracAligned');
    disp(S);
end