function reportMultiFaultContribution(faultResults, subsectionData, invData, bestRates)
% reportMultiFaultContribution - Calculates the absolute and percentage 
% contribution of single-fault vs multi-fault ruptures in terms of both 
% event rates (frequency) and seismic moment rates, and plots a stacked MFD chart.

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
    
    % Pre-allocate vectors for rupture properties
    numActive = length(activeRupIdx);
    rupMw = zeros(numActive, 1);
    isMultiFault = false(numActive, 1);
    
    logicalGsr = invData.Gsr > 0;
    
    counter = 0;
    for r = activeRupIdx'
        counter = counter + 1;
        
        % Identify involved subsections and map to parent faults
        involvedSubs = find(logicalGsr(r, :));
        involvedFaults = unique(subToFaultMap(involvedSubs));
        involvedFaults(involvedFaults == 0) = []; % Clean padding
        
        % Retrieve Rupture Magnitude (Mw)
        if isfield(invData, 'Mw_r')
            rupMw(counter) = invData.Mw_r(r);
        elseif isfield(invData, 'Mw')
            rupMw(counter) = invData.Mw(r);
        else
            % Fallback: compute Mw from Moment if Mw vector is not directly in invData
            % M0 = mu * Area * Slip
            if isfield(invData, 'Area_r')
                rupArea_m2 = invData.Area_r(r); 
            else
                rupArea_m2 = sum(subsectionData.Area_km2(involvedSubs)) * 1e6;
            end
            M0 = mu * rupArea_m2 * invData.Dr(r);
            rupMw(counter) = (log10(M0) - 9.05) / 1.5;
        end
        
        % Calculate Seismic Moment Rate
        if isfield(invData, 'Area_r')
            rupArea_m2 = invData.Area_r(r); 
        else
            rupArea_m2 = sum(subsectionData.Area_km2(involvedSubs)) * 1e6;
        end
        
        rupMoRate = bestRates.mean(r) * invData.Dr(r) * rupArea_m2 * mu;
        rupRate = bestRates.mean(r);
        
        % Classify based on how many parent faults are involved
        if length(involvedFaults) > 1
            % Multi-Fault Rupture
            isMultiFault(counter) = true;
            totalRateMulti = totalRateMulti + rupRate;
            totalMoRateMulti = totalMoRateMulti + rupMoRate;
        else
            % Single-Fault Rupture
            isMultiFault(counter) = false;
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
    fprintf('📊 SEISMIC BUDGETING: SINGLE VS MULTI-FAULT EVENTS\n');
    fprintf('=========================================================\n');
    fprintf('Total Active Ruptures in Solution: %d\n\n', numActive);
    
    fprintf('📈 BY EVENT FREQUENCY RATE (events/yr):\n');
    fprintf('  ▪️ Single-Fault Ruptures:  %.5f/yr  (%.1f%%)\n', totalRateSingle, pctRateSingle);
    fprintf('  ▪️ Multi-Fault Ruptures:   %.5f/yr  (%.1f%%)\n', totalRateMulti, pctRateMulti);
    fprintf('  Total System Event Rate:   %.5f/yr\n\n', grandTotalRate);
    
    fprintf('⚡ BY SEISMIC MOMENT RATE RELEASE (N·m/yr):\n');
    fprintf('  ▪️ Single-Fault Ruptures:  %.2e N·m/yr (%.1f%%)\n', totalMoRateSingle, pctMoRateSingle);
    fprintf('  ▪️ Multi-Fault Ruptures:   %.2e N·m/yr (%.1f%%)\n', totalMoRateMulti, pctMoRateMulti);
    fprintf('  Total System Moment Rate:  %.2e N·m/yr\n', grandTotalMoRate);
    fprintf('=========================================================\n\n');
    
    % =====================================================================
    % PLOT: STACKED HISTOGRAM OF RUPTURE COUNTS BY MAGNITUDE (Mw)
    % =====================================================================
    mwBinWidth = 0.1; % Magnitude bin size (e.g., 0.1 or 0.2)
    minMw = floor(min(rupMw) / mwBinWidth) * mwBinWidth;
    maxMw = ceil(max(rupMw) / mwBinWidth) * mwBinWidth;
    mwEdges = minMw:mwBinWidth:maxMw;
    
    % Bin active ruptures by magnitude and rupture type
    singleCounts = histcounts(rupMw(~isMultiFault), mwEdges);
    multiCounts  = histcounts(rupMw(isMultiFault), mwEdges);
    
    % Bin centers for x-axis plotting
    mwCenters = mwEdges(1:end-1) + mwBinWidth/2;
    
    figure('Name', 'Rupture Distribution by Magnitude', 'Color', 'w');
    
    % Create Stacked Bar Chart: [Single-Fault, Multi-Fault]
    hBar = bar(mwCenters, [singleCounts', multiCounts'], 1.0, 'stacked');
    
    % Set colors: Blue for Single-Fault, Red for Multi-Fault
    hBar(1).FaceColor = [0.2 0.4 0.8];
    hBar(2).FaceColor = [0.85 0.25 0.2];
    
    grid on;
    xlabel('Magnitude (M_w)', 'Interpreter', 'tex');
    ylabel('Number of Ruptures');
    title('Single-Fault vs Multi-Fault ruptures by magnitude bin');
    legend({'Single-Fault Ruptures', 'Multi-Fault Ruptures'}, 'Location', 'northeast');
    
    % Add total text on top of high bars for extra clarity
    ylimits = ylim;
    ylim([0 ylimits(2) * 1.05]);
end
