function reportMultiFaultContribution(faultResults, subsectionData, invData, bestRates)
% reportMultiFaultContribution - Calculates the absolute and percentage 
% contribution of single-fault vs multi-fault ruptures in terms of both 
% event rates (frequency) and seismic moment rates.

    mu = 3e10; % Rigidity / Shear Modulus (Pa)
    numFaults = length(faultResults);
    activeRupIdx = find(bestRates.mean > 0);
    
    % Map Global Subsection ID -> Parent Fault ID
    subToFaultMap = zeros(height(subsectionData), 1);
    for nf = 1:numFaults
        subIdx = faultResults(nf).GlobalSubIndices;
        subToFaultMap(subIdx) = nf;
    end
    
    % Initialize counters for Rates and Moment Rates
    totalRateSingle = 0;
    totalRateMulti  = 0;
    totalMoRateSingle = 0;
    totalMoRateMulti  = 0;
    
    logicalGsr = invData.Gsr > 0;
    
    for r = activeRupIdx'
        % Identify involved subsections and map to parent faults
        involvedSubs = find(logicalGsr(r, :));
        involvedFaults = unique(subToFaultMap(involvedSubs));
        involvedFaults(involvedFaults == 0) = []; % Clean padding
        
        % Calculate the seismic moment rate of this specific rupture:
        % MoRate = rate * slip * area * mu
        % (invData.Area_r potrebbe essere disponibile, altrimenti lo ricaviamo dalle sottosezioni)
        % Se invData ha già l'area della rottura usiamo quella, altrimenti sommiamo l'area delle sottosezioni coinvolte
        if isfield(invData, 'Area_r')
            rupArea_m2 = invData.Area_r(r); 
        else
            rupArea_m2 = sum(subsectionData.Area_km2(involvedSubs)) * 1e6; % Convert km2 to m2
        end
        
        rupMoRate = bestRates.mean(r) * invData.Dr(r) * rupArea_m2 * mu;
        rupRate = bestRates.mean(r);
        
        % Classify based on how many parent faults are involved
        if length(involvedFaults) > 1
            % Multi-Fault Rupture
            totalRateMulti = totalRateMulti + rupRate;
            totalMoRateMulti = totalMoRateMulti + rupMoRate;
        else
            % Single-Fault Rupture
            totalRateSingle = totalRateSingle + rupRate;
            totalMoRateSingle = totalMoRateSingle + rupMoRate;
        end
    end
    
    % Sum totals
    grandTotalRate = totalRateSingle + totalRateMulti;
    grandTotalMoRate = totalMoRateSingle + totalMoRateMulti;
    
    % Calculate Percentages
    pctRateSingle = (totalRateSingle / grandTotalRate) * 100;
    pctRateMulti  = (totalRateMulti / grandTotalRate) * 100;
    
    pctMoRateSingle = (totalMoRateSingle / grandTotalMoRate) * 100;
    pctMoRateMulti  = (totalMoRateMulti / grandTotalMoRate) * 100;
    
    % Print clean report to command window
    fprintf('\n=========================================================\n');
    fprintf('📊 SEISMIC BUDGET BUDGETING: SINGLE VS MULTI-FAULT EVENTS\n');
    fprintf('=========================================================\n');
    fprintf('Total Active Ruptures in Solution: %d\n\n', length(activeRupIdx));
    
    fprintf('📈 BY EVENT FREQUENCY RATE (events/yr):\n');
    fprintf('  ▪️ Single-Fault Ruptures:  %.5f/yr  (%.1f%%)\n', totalRateSingle, pctRateSingle);
    fprintf('  ▪️ Multi-Fault Ruptures:   %.5f/yr  (%.1f%%)\n', totalRateMulti, pctRateMulti);
    fprintf('  Total System Event Rate:   %.5f/yr\n\n', grandTotalRate);
    
    fprintf('⚡ BY SEISMIC MOMENT RATE RELEASE (N·m/yr):\n');
    fprintf('  ▪️ Single-Fault Ruptures:  %.2e N·m/yr (%.1f%%)\n', totalMoRateSingle, pctMoRateSingle);
    fprintf('  ▪️ Multi-Fault Ruptures:   %.2e N·m/yr (%.1f%%)\n', totalMoRateMulti, pctMoRateMulti);
    fprintf('  Total System Moment Rate:  %.2e N·m/yr\n', grandTotalMoRate);
    fprintf('=========================================================\n\n');
end