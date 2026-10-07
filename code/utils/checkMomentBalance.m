function checkMomentBalance(bestRates, invData, subsectionData, targets, params)
% CHECKMOMENTBALANCE - Compares tectonic moment range vs. model release per fault.
    fprintf('⚖️ Calculating Moment Rate Budget Allocation...\n');
    
    mu = params.mu; 
    faultNames = unique(subsectionData.FaultName, 'stable');
    numFaults = numel(faultNames);
    
    budgetMin   = zeros(numFaults, 1);
    budgetMax   = zeros(numFaults, 1);
    budgetMean  = zeros(numFaults, 1);
    modelTotal  = zeros(numFaults, 1);
    
    % Logical structure of Gsr for efficient mapping
    logicalGsr = invData.Gsr > 0; 
    activeRupIdx = find(bestRates.mean > 1e-12);
    
    for i = 1:numFaults
        fName = faultNames{i};
        subIdx = find(strcmp(subsectionData.FaultName, fName));
        
        % --- 1. Geological Budget Range (Target Min, Max, Mean) ---
        areas_m2 = subsectionData.Area_km2(subIdx) * 1e6;
        
        budgetMin(i)  = sum(areas_m2 .* targets.minSR(subIdx)) * mu;
        budgetMax(i)  = sum(areas_m2 .* targets.maxSR(subIdx)) * mu;
        budgetMean(i) = (budgetMin(i) + budgetMax(i)) / 2;
        
        % --- 2. Model Output (Inversion Result) ---
        faultMoment = 0;
        for r = activeRupIdx'
            subsInRupture = find(logicalGsr(r, :));
            commonSubs = intersect(subIdx, subsInRupture);
            
            if ~isempty(commonSubs)
                areaOnFault_m2 = sum(subsectionData.Area_km2(commonSubs)) * 1e6;
                slip = invData.Dr(r);
                faultMoment = faultMoment + bestRates.mean(r) * slip * areaOnFault_m2 * mu;
            end
        end
        modelTotal(i) = faultMoment;
    end
    
    % --- 3. Visualization ---
    figure('Name', 'Moment Rate Allocation', 'Color', 'w');
    
    x = 1:numFaults;
    
    % A. Single bar for Inversion Result
    b = bar(x, modelTotal, 0.5, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', ...
            'DisplayName', 'Inversion (Mean)');
    hold on;
    
    % B. Errorbar representing [Min, Max] Geological Range
    % Lower error = budgetMean - budgetMin, Upper error = budgetMax - budgetMean
    errLow  = budgetMean - budgetMin;
    errHigh = budgetMax - budgetMean;
    
    hErr = errorbar(x, budgetMean, errLow, errHigh, 'k.', 'LineWidth', 1.2, ...
                    'CapSize', 5, 'DisplayName', 'Tectonic Moment Rate Range [Min, Max]');
                
    % C. Square marker representing Geological Mean Target
    hSq = plot(x, budgetMean, 's', 'MarkerSize', 6, ...
               'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.2 0.4 0.6], ...
               'LineWidth', 1.0, 'DisplayName', 'Tectonic Moment Rate (Mean)');
    
    % Formatting
    set(gca, 'XTick', 1:numFaults, 'XTickLabel', faultNames, 'XTickLabelRotation', 45);
    ylabel('Moment Rate [N-m/yr]', 'FontSize', 11);
    title('Seismic Moment Rate: Tectonic vs. Inversion', 'FontSize', 12);
    legend([b, hSq, hErr], 'Location', 'northeast');
    grid on; box on;
    xlim([0.5, numFaults + 0.5]);
    
    % --- 4. Final Statistics ---
    totalBudget = sum(budgetMean);
    totalModel  = sum(modelTotal);
    
    fprintf('\n--- Global Moment Report ---\n');
    fprintf('Total Tectonic Budget (Mean): %.2e Nm/yr\n', totalBudget);
    fprintf('Total Model Release:            %.2e Nm/yr\n', totalModel);
    fprintf('📊 Global Balance: Model covers %.2f%% of the Mean Tectonic Budget.\n', (totalModel/totalBudget)*100);
    
    if (totalModel/totalBudget) < 0.8
        fprintf('⚠️ Warning: Model could under-predicting the total moment budget.\n');
    elseif (totalModel/totalBudget) > 1.2
        fprintf('⚠️ Warning: Model could over-predicting the total moment budget.\n');
    end
end
