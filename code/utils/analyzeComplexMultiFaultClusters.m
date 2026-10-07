function complexRuptureStats = analyzeComplexMultiFaultClusters(faultResults, subsectionData, invData, bestRates)
% analyzeComplexMultiFaultClusters - Identifies and ranks unique clusters 
% of 3 or more parent faults breaking together in a single seismic event,
% and visualizes the distribution of faults involved per multi-fault rupture.

    numFaults = length(faultResults);
    activeRupIdx = find(bestRates.mean > 0);
    faultNames = {faultResults.Name};
    
    % Map Subsection Global ID -> Parent Fault ID
    subToFaultMap = zeros(height(subsectionData), 1);
    for nf = 1:numFaults
        subIdx = faultResults(nf).GlobalSubIndices;
        subToFaultMap(subIdx) = nf;
    end
    
    % Container for 3+ fault clusters
    clusterMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Vector to store the number of involved parent faults per active rupture
    numActive = length(activeRupIdx);
    faultsPerRupture = zeros(numActive, 1);
    
    logicalGsr = invData.Gsr > 0;
    fprintf('🔍 Scanning for high-order multi-fault ruptures (3+ faults)...\n');
    
    counter = 0;
    for r = activeRupIdx'
        counter = counter + 1;
        involvedSubs = find(logicalGsr(r, :));
        involvedFaults = unique(subToFaultMap(involvedSubs));
        involvedFaults(involvedFaults == 0) = []; % Remove padding
        
        % Store count of parent faults involved in this active rupture
        faultsPerRupture(counter) = length(involvedFaults);
        
        % Filter only strict clusters of 3 or more faults
        if length(involvedFaults) >= 3
            % Sort IDs to ensure (A,B,C) is the same as (C,A,B)
            sortedFaults = sort(involvedFaults);
            
            % Create a unique string key (e.g., '1_4_12')
            key = strjoin(arrayfun(@num2str, sortedFaults, 'UniformOutput', false), '_');
            
            if clusterMap.isKey(key)
                data = clusterMap(key);
                data.RuptureCount = data.RuptureCount + 1;
                data.CombinedRate = data.CombinedRate + bestRates.mean(r);
                clusterMap(key) = data;
            else
                data = struct();
                data.FaultIDs = sortedFaults;
                data.RuptureCount = 1;
                data.CombinedRate = bestRates.mean(r);
                clusterMap(key) = data;
            end
        end
    end
    
    % Convert the Map to a readable Table
    keys = clusterMap.keys();
    if isempty(keys)
        fprintf('ℹ️ No ruptures involving 3 or more faults were activated by the inversion.\n');
        complexRuptureStats = table();
    else
        clusterList = [];
        for k = 1:length(keys)
            data = clusterMap(keys{k});
            
            % Resolve names of all faults involved
            namesInCluster = faultNames(data.FaultIDs);
            clusterString = strjoin(namesInCluster, ' + ');
            
            clusterList = [clusterList; struct(...
                'FaultCluster', clusterString, ...
                'NumberOfFaults', length(data.FaultIDs), ...
                'UniqueRuptures', data.RuptureCount, ...
                'CombinedRate', data.CombinedRate, ...
                'MeanReturnPeriod', 1 / data.CombinedRate)];
        end
        
        % Format and Sort Output Table
        T = struct2table(clusterList);
        T = sortrows(T, 'CombinedRate', 'descend');
        
        fprintf('\n🌋 TOP COMPLEX MULTI-FAULT CLUSTERS (Ranked by Rate):\n');
        disp(T(1:min(10, height(T)), :));
        
        complexRuptureStats = T;
        
        % =====================================================================
        % PLOT 1: TOP 10 COMPLEX MULTI-FAULT CLUSTERS (3+ FAULTS)
        % =====================================================================
        figure('Name', '3+ Multifault Rupture Rates', 'Color', 'w');
        numBars = min(10, height(T));
        barh(categorical(T.FaultCluster(numBars:-1:1)), T.CombinedRate(numBars:-1:1), ...
             'FaceColor', [0.6350 0.0780 0.1840]);
        xlabel('Event Rate (events/yr)');
        % title('Top Complex Seismic Cascades (3+ Faults)');
        grid on;
    end
    
    % =====================================================================
    % PLOT 2: NUMBER OF FAULTS INVOLVED PER MULTIFAULT RUPTURE (>= 2 FAULTS)
    % =====================================================================
    multiFaultRups = faultsPerRupture(faultsPerRupture >= 2);
    
    if ~isempty(multiFaultRups)
        % Get unique levels of fault counts and their occurrence frequencies
        [uniqueFaultCounts, ~, binIdx] = unique(multiFaultRups);
        ruptureCountsPerLevel = accumarray(binIdx, 1);
        
        figure('Name', 'Faults Involved per Multi-Fault Rupture', 'Color', 'w');
        
        % Plot bar chart
        bar(categorical(uniqueFaultCounts), ruptureCountsPerLevel, ...
            'FaceColor', [0.85 0.325 0.098], 'EdgeColor', 'k');
        
        grid on;
        xlabel('Number of faults involved in rupture');
        ylabel('Number of multifault ruptures');
        title('How many faults are involved per multifault rupture?', ...
              'FontSize', 12, 'FontWeight', 'bold');
        
        % Add exact value above each bar
        ylimits = ylim;
        ylim([0 ylimits(2) * 1.08]); % Leave headroom for labels
        
        for b = 1:length(uniqueFaultCounts)
            text(b, ruptureCountsPerLevel(b) + ylimits(2)*0.015, ...
                 num2str(ruptureCountsPerLevel(b)), ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', ...
                 'FontSize', 9, 'FontWeight', 'normal');
        end
        
        xtickangle(45); % Rotate x-axis labels for optimal readability
    end
end
