function complexRuptureStats = analyzeComplexMultiFaultClusters(faultResults, subsectionData, invData, bestRates)
% analyzeComplexMultiFaultClusters - Identifies and ranks unique clusters 
% of 3 or more parent faults breaking together in a single seismic event.

    numFaults = length(faultResults);
    activeRupIdx = find(bestRates.mean > 0);
    faultNames = {faultResults.Name};
    
    % Map Subsection Global ID -> Parent Fault ID
    subToFaultMap = zeros(height(subsectionData), 1);
    for nf = 1:numFaults
        subIdx = faultResults(nf).GlobalSubIndices;
        subToFaultMap(subIdx) = nf;
    end
    
    % Container for our clusters
    % We will use a map where the key is the sorted combination of fault IDs
    clusterMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    logicalGsr = invData.Gsr > 0;
    fprintf('🔍 Scanning for high-order multi-fault ruptures (3+ faults)...\n');
    
    for r = activeRupIdx'
        involvedSubs = find(logicalGsr(r, :));
        involvedFaults = unique(subToFaultMap(involvedSubs));
        involvedFaults(involvedFaults == 0) = [];
        
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
        return;
    end
    
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
    
    % --- VISUALIZATION ---
    figure('Name', 'NextQuake - 3+ Multi-Fault Cluster Rates', 'Color', 'w');
    numBars = min(10, height(T));
    barh(categorical(T.FaultCluster(numBars:-1:1)), T.CombinedRate(numBars:-1:1), 'FaceColor', [0.6350 0.0780 0.1840]);
    xlabel('Combined Event Rate (events/yr)');
    title('Top Complex Seismic Cascades (3+ Faults)');
    grid on;
end