function checkMomentBalance(bestRates, invData, subsectionData, targets, params)
% CHECKMOMENTBALANCE - Compares the tectonic moment budget vs. the model release.
% This version uses the mean of the box-constraints (minSR/maxSR) as the target.

    fprintf('⚖️ Calculating Moment Rate Balance (using Gsr matrix)...\n');
    
    mu = params.mu; 
    faultNames = unique(subsectionData.FaultName, 'stable');
    numFaults = numel(faultNames);
    
    budgetTotal = zeros(numFaults, 1);
    modelTotal = zeros(numFaults, 1);
    
    % Use the logical structure of Gsr for efficient subsection-rupture mapping
    logicalGsr = invData.Gsr > 0; 
    % Optimize: only loop through ruptures with a non-zero rate
    activeRupIdx = find(bestRates.mean > 1e-12);

    for i = 1:numFaults
        fName = faultNames{i};
        % Get global indices of subsections belonging to the current fault
        subIdx = find(strcmp(subsectionData.FaultName, fName));
        
        % --- 1. Geological Budget (Target) ---
        % Since we use Box-Constraints, we define the target as the mid-point 
        % of the allowed Slip Rate range [minSR, maxSR].
        meanTargetSR = (targets.minSR(subIdx) + targets.maxSR(subIdx)) / 2;
        
        % Budget [Nm/yr] = Area [m^2] * SlipRate [m/yr] * mu [Pa]
        % Note: Area_km2 is converted to m2 by multiplying by 1e6
        budgetTotal(i) = sum(subsectionData.Area_km2(subIdx) * 1e6 .* meanTargetSR) * mu;
        
        % --- 2. Model Output (Inversion Result) ---
        faultMoment = 0;
        
        for r = activeRupIdx'
            % Find subsections involved in the current rupture 'r'
            subsInRupture = find(logicalGsr(r, :));
            
            % Identify subsections that belong BOTH to this fault and this rupture
            commonSubs = intersect(subIdx, subsInRupture);
            
            if ~isempty(commonSubs)
                % Calculate the area of the portion of the fault that slipped
                areaOnFault_m2 = sum(subsectionData.Area_km2(commonSubs)) * 1e6;
                % Average slip of the rupture
                slip = invData.Dr(r);
                
                % Increment moment rate for this fault: Rate * Slip * Area * mu
                faultMoment = faultMoment + bestRates.mean(r) * slip * areaOnFault_m2 * mu;
            end
        end
        modelTotal(i) = faultMoment;
    end
    
    % --- 3. Visualization ---
    figure('Name', 'Moment Rate Balance', 'Color', 'w', 'Position', [100 100 900 500]);
    b = bar([budgetTotal, modelTotal], 'EdgeColor', 'none');
    b(1).FaceColor = [0.2 0.4 0.6]; % Steel Blue for Geological Budget
    b(2).FaceColor = [0.8 0.3 0.3]; % Terracotta for Inversion Result
    
    set(gca, 'XTick', 1:numFaults, 'XTickLabel', faultNames, 'XTickLabelRotation', 45);
    ylabel('Moment Rate [N-m/yr]');
    legend({'Geological Budget (Mean)', 'Inversion Result'}, 'Location', 'northeast');
    grid on;
    title('Moment Balance: Tectonic Budget vs. Model Release');
    
    % --- 4. Final Statistics ---
    totalBudget = sum(budgetTotal);
    totalModel = sum(modelTotal);
    
    fprintf('\n--- Global Moment Report ---\n');
    fprintf('Total Budget: %.2e Nm/yr\n', totalBudget);
    fprintf('Total Model:  %.2e Nm/yr\n', totalModel);
    fprintf('📊 Global Balance: Model covers %.2f%% of the Geological Budget.\n', (totalModel/totalBudget)*100);
    
    if (totalModel/totalBudget) < 0.8
        fprintf('⚠️ Warning: Model is significantly under-predicting the moment. Check weights constraints.\n');
    elseif (totalModel/totalBudget) > 1.2
        fprintf('⚠️ Warning: Model is over-predicting the moment. Check weights constraints.\n');
    end
end