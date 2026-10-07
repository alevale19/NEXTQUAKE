function faultCoRuptureStats = analyzeFaultCoRupturePairs(faultResults, subsectionData, invData, bestRates)
% analyzeFaultCoRupturePairs - Analyzes and visualizes how frequently 
% parent faults rupture together in multi-fault inversion solutions.

    numFaults = length(faultResults);
    activeRupIdx = find(bestRates.mean > 0);
    
    % 1. Create a mapping vector: Subsection Global ID -> Parent Fault ID
    subToFaultMap = zeros(height(subsectionData), 1);
    for nf = 1:numFaults
        subIdx = faultResults(nf).GlobalSubIndices;
        subToFaultMap(subIdx) = nf;
    end
    
    % 2. Initialize Matrices for Frequencies and Rates
    % Matrix of shared event rates (events/yr)
    coRuptureRates = zeros(numFaults, numFaults);
    % Matrix of actual counts of active ruptures shared
    coRuptureCounts = zeros(numFaults, numFaults);
    
    fprintf('🧮 Analyzing multi-fault co-rupture statistics...\n');
    
    % 3. Loop over each active rupture in the solution
    logicalGsr = invData.Gsr > 0;
    for r = activeRupIdx'
        % Find which global subsections are involved in this specific rupture
        involvedSubs = find(logicalGsr(r, :));
        
        % Map these subsections to their Parent Fault IDs
        involvedFaults = unique(subToFaultMap(involvedSubs));
        involvedFaults(involvedFaults == 0) = []; % Remove any unmapped padding
        
        % If the rupture involves more than 1 parent fault, it's a multi-fault event!
        if length(involvedFaults) > 1
            % Permute all pairs of faults involved in this rupture
            for i = 1:length(involvedFaults)
                for j = 1:length(involvedFaults)
                    if i ~= j % Avoid self-pairing on diagonal for co-rupture analysis
                        f1 = involvedFaults(i);
                        f2 = involvedFaults(j);
                        
                        % Accumulate rates and counts
                        coRuptureRates(f1, f2) = coRuptureRates(f1, f2) + bestRates.mean(r);
                        coRuptureCounts(f1, f2) = coRuptureCounts(f1, f2) + 1;
                    end
                end
            end
        end
    end
    
    % Since the matrices are symmetric, let's extract the upper triangle for ranking
    pairList = [];
    faultNames = {faultResults.Name};
    
    for i = 1:numFaults
        for j = i+1:numFaults
            if coRuptureCounts(i, j) > 0
                pairList = [pairList; struct(...
                    'Fault1', faultNames{i}, ...
                    'Fault2', faultNames{j}, ...
                    'RuptureCount', coRuptureCounts(i, j), ...
                    'CombinedRate', coRuptureRates(i, j), ...
                    'MeanReturnPeriod', 1 / coRuptureRates(i, j))];
            end
        end
    end
    
    % 4. Sort and Print the Top Co-Rupturing Pairs by Rate
    if ~isempty(pairList)
        T = struct2table(pairList);
        T = sortrows(T, 'CombinedRate', 'descend');
        
        fprintf('\n🏆 TOP MULTI-FAULT CO-RUPTURE PAIRS (Ranked by Rate):\n');
        disp(T(1:min(15, height(T)), :));
        faultCoRuptureStats.Table = T;
    else
        fprintf('ℹ️ No multi-fault co-ruptures found active in the current inversion solution.\n');
        faultCoRuptureStats.Table = table();
    end
    
    % 5. PLOT THE CO-RUPTURE MATRICES
    figure('Name', 'Fault Section Co-Rupture Matrix', 'Color', 'w');
    
    % Subplot 1: Shared Rates (Log10)
    subplot(1, 2, 1);
    logRatesPlot = log10(coRuptureRates);
    logRatesPlot(coRuptureRates == 0) = NaN; % Hide zeros
    imagesc(logRatesPlot);
    colormap(gca, jet); cb1 = colorbar; 
    ylabel(cb1, 'Co-Rupture Rate [Log_{10}(events/yr)]');
    set(gca, 'XTick', 1:numFaults, 'XTickLabel', faultNames, 'XTickLabelRotation', 45, ...
             'YTick', 1:numFaults, 'YTickLabel', faultNames, 'FontSize', 8, 'TickLabelInterpreter', 'none');
    title('Shared Annual Event Rates');
    axis square; grid on;
    
    % Subplot 2: Actual Rupture Counts
    subplot(1, 2, 2);
    imagesc(coRuptureCounts);
    colormap(gca, flipud(hot)); cb2 = colorbar; 
    ylabel(cb2, 'Number of Active Unique Ruptures');
    set(gca, 'XTick', 1:numFaults, 'XTickLabel', faultNames, 'XTickLabelRotation', 45, ...
             'YTick', 1:numFaults, 'YTickLabel', faultNames, 'FontSize', 8, 'TickLabelInterpreter', 'none');
    title('Shared Unique Active Ruptures (Count)');
    axis square; grid on;
    
    % Save data matrices into output struct
    faultCoRuptureStats.RatesMatrix = coRuptureRates;
    faultCoRuptureStats.CountsMatrix = coRuptureCounts;
end
