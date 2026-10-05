function [bestRates, energyHistory] = runSimulatedAnnealing(invData, targets, params, faultResults)
% RUNSIMULATEDANNEALING - Classic Parallel Ensemble with Seismotectonic Delta Calibration
% Solves for rupture rates by minimizing total energy using a static deltaScale 
% anchored to the regional minimum plausible geological slip rate.

    % --- 1. HYPERPARAMETERS & SEED SETUP ---
    numIter = params.numIter;
    numRuns = params.runSA;
    T_init  = 10.0;
    T_final = 0.001;
    cooling_rate = (T_final / T_init)^(1 / numIter); 
    
    numRup = invData.numRuptures;
    numFaults = length(faultResults);
    
    % --- SEISMOTECTONIC STATIC DELTASCALE CALIBRATION ---
    % Find positive slip rates (exclude unconstrained or zero-rate sections)
    validSR = targets.maxSR(targets.maxSR > 0);
    if ~isempty(validSR)
        % Reference lower-bound slip rate (10th percentile to avoid unrepresentative outliers)
        referenceSR = prctile(validSR, 10); 
        % Convert to standard return period rate assuming a ~1m typical event slip displacement
        % Rate = 1 / ReturnPeriod = SlipRate / 1m
        staticDeltaScale = referenceSR; 
    else
        staticDeltaScale = 1e-4; % Safe fallback for regional crustal faults
    end
    
    % Filtering valid ruptures (exclude segmentation violations)
    validRupIdx = setdiff(1:numRup, invData.badRupIdx);
    numValid = numel(validRupIdx);
    
    % =========================================================================
    % ⚠️ INTELLIGENT WARNING: INITIAL ITERATION SUFFICIENCY CHECK
    % =========================================================================
    % 1. Calculate the average number of samplings each valid rupture receives
    samplingFactorK = numIter / numValid;
    
    % Engineering critical thresholds
    minRecommendedK = params.minRecommendedK;
    
    if samplingFactorK < minRecommendedK
        warning('⚠️ INSUFFICIENT ITERATIONS DETECTED relative to the scale of the problem!');
        fprintf('   ---------------------------------------------------------------------------\n');
        fprintf('   👉 You have set %d iterations for %d valid ruptures.\n', numIter, numValid);
        fprintf('   👉 Each rupture will be sampled an average of %.1f times.\n', samplingFactorK);
        fprintf('   🔔 SUGGESTION: Increase the number of iterations to at least %d for this configuration,\n', numValid * 1000);
        fprintf('      or proceed knowing that the results might be underfitted.\n');
        fprintf('   ---------------------------------------------------------------------------\n\n');
        
        % 3-second courtesy pause to give the user time to read the warning before parfor starts
        pause(3); 
    end
    % =========================================================================
    % Pre-allocate ensemble variables for parfor
    ensembleRates = zeros(numRup, numRuns); 
    energyHistory = zeros(numIter, numRuns); 
    faultSlipEnsemble = zeros(numFaults, numRuns);
    
    fprintf('\n🚀 Starting Parallel SA Ensemble (%d Runs) with perturbation of: %.2e...\n', numRuns, staticDeltaScale);

    % --- 2. PARALLEL ENSEMBLE LOOP ---
    parfor n = 1:numRuns
        fprintf('▶️ Worker: starting Run %d...\n', n);
        
        % Reproducibility anchor per worker
        if isfield(params, 'seed')
            rng(params.seed + n);
        end
        
        % Initialization
        localRates = zeros(numRup, 1);
        [currentEnergy, ~] = calculateEnergy(localRates, invData, targets, params);
        
        bestRatesInRun = localRates;
        bestEnergyInRun = currentEnergy;
        T = T_init;
        
        % --- Main Optimization Loop ---
        for i = 1:numIter
            % 1. Proposal
            idx = validRupIdx(randi(numValid));
            oldRate = localRates(idx);
            
            % 2. Perturbation: Static physical step size
            delta = (rand() - 0.5) * staticDeltaScale;
            newRate = max(0, oldRate + delta); 
            
            % 3. Energy Calculation
            tempRates = localRates;
            tempRates(idx) = newRate;
            [newEnergy, ~] = calculateEnergy(tempRates, invData, targets, params);
            
            % 4. Metropolis Acceptance Criterion
            deltaE = newEnergy - currentEnergy;
            if deltaE < 0 || rand() < exp(-deltaE / T)
                localRates(idx) = newRate;
                currentEnergy = newEnergy;
                
                if currentEnergy < bestEnergyInRun
                    bestEnergyInRun = currentEnergy;
                    bestRatesInRun = localRates;
                end
            end
            
            % 5. Cooling Schedule
            T = T * cooling_rate;
            energyHistory(i, n) = currentEnergy; 
        end
        
        ensembleRates(:, n) = bestRatesInRun; 
        
        % --- BRIDGE: Final Fault-Level Slip Rate Calculation ---
        currentSR_sub = (invData.Gsr' * (bestRatesInRun .* invData.Dr));
        localFaultSlip = zeros(numFaults, 1);
        for fIdx = 1:numFaults
            subIdx = faultResults(fIdx).GlobalSubIndices;
            if ~isempty(subIdx)
                localFaultSlip(fIdx) = mean(currentSR_sub(subIdx));
            end
        end
        faultSlipEnsemble(:, n) = localFaultSlip; 
        
        fprintf('✅ Worker: Run %d finished. Final Energy: %.2e\n', n, bestEnergyInRun);
    end
    
    % --- 3. EPISTEMIC ENSEMBLE ANALYSIS ---
    [~, absoluteBestRunIdx] = min(energyHistory(end, :));
    
    bestRates = struct();

    bestRates.matrix     = ensembleRates;
    
    bestRates.mean       = mean(ensembleRates, 2);
    bestRates.p5         = prctile(ensembleRates, 5, 2);
    bestRates.p16        = prctile(ensembleRates, 16, 2);
    bestRates.p50        = prctile(ensembleRates, 50, 2);
    bestRates.p84        = prctile(ensembleRates, 84, 2);
    bestRates.p95        = prctile(ensembleRates, 95, 2);
    bestRates.bestOfBest = ensembleRates(:, absoluteBestRunIdx);
    
    if isfield(params, 'seed')
        bestRates.runSeeds = params.seed + (1:numRuns)';
        bestRates.bestOfBestSeed = params.seed + absoluteBestRunIdx;
    end
    
    % % --- 4. STABILITY REPORT ---
    % fprintf('\n--- Final Stability Report (Ensemble Analysis) ---\n');
    % mean_rup = bestRates.mean;
    % active_rup = mean_rup > 1e-10;
    % if any(active_rup)
    %     std_rup = std(ensembleRates(active_rup, :), 0, 2);
    %     cv_rup = mean(std_rup ./ mean_rup(active_rup)) * 100;
    %     %fprintf('📊 Rupture-level CV: %.2f%%\n', cv_rup);
    % end
    % 
    % mean_fault = mean(faultSlipEnsemble, 2);
    % active_fault = mean_fault > 1e-6;
    % if any(active_fault)
    %     std_fault = std(faultSlipEnsemble(active_fault, :), 0, 2);
    %     cv_fault = mean(std_fault ./ mean_fault(active_fault)) * 100;
    %     %fprintf('⭐ Fault-level Slip Rate CV: %.2f%%\n', cv_fault);
    % end
    fprintf('🏆 Global Minimum Solution Found in Run #%d.\n', absoluteBestRunIdx);

    % 1. Compute dimensionless Coefficient of Variation (CV = std / mean)
mean_rup_all = mean(ensembleRates, 2);
std_rup_all  = std(ensembleRates, 0, 2);
active_rup   = mean_rup_all > 1e-6; % Filter out inactive/negligible ruptures
cv_rup_vec   = std_rup_all(active_rup) ./ mean_rup_all(active_rup); % Dimensionless CV

% Compute fault-level slip rates across all ensemble runs
all_runs_SR  = invData.Gsr' * (ensembleRates .* invData.Dr);
mean_fault   = mean(all_runs_SR, 2);
std_fault    = std(all_runs_SR, 0, 2);
active_fault = mean_fault > 1e-6;
cv_fault_vec = std_fault(active_fault) ./ mean_fault(active_fault); % Dimensionless CV

% 2. Plot Ensemble CV Histograms
figure('Name', 'Ensemble CV Histograms', 'Position', [100 100 1000 420]);

% Panel A: Rupture-Level Rate Variability
subplot(1, 2, 1);
histogram(cv_rup_vec, 30, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w');
grid on; box on;
title('a) Rupture-Level Rate Variability', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Coefficient of Variation (CV)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Number of Ruptures', 'FontSize', 10, 'FontWeight', 'bold');
xline(median(cv_rup_vec), 'r--', sprintf(' Median: %.2f', median(cv_rup_vec)), ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 10);

% Panel B: Subsection-Level Slip Rate Variability
subplot(1, 2, 2);
histogram(cv_fault_vec, 15, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'w');
grid on; box on;
title('b) Subsection-Level Slip Rate Variability', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Coefficient of Variation (CV)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Number of Subsections', 'FontSize', 10, 'FontWeight', 'bold');
xline(median(cv_fault_vec), 'r--', sprintf(' Median: %.2f', median(cv_fault_vec)), ...
      'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 10);

sgtitle('Ensemble Solution Variability Across Independent Runs', ...
        'FontSize', 12, 'FontWeight', 'bold');

% 1. Project all ensemble runs to subsection slip rates [numSubsections x numRuns]
all_runs_SR = invData.Gsr' * (ensembleRates .* invData.Dr);

% 2. Calculate ensemble statistics directly on subsection slip rates
meanSR = mean(all_runs_SR, 2);
p5_SR  = prctile(all_runs_SR, 5, 2);
p95_SR = prctile(all_runs_SR, 95, 2);

numSub = length(meanSR);
x_idx  = 1:numSub;

% Fig. 2 - Slip Rate Match with 5th-95th Percentile Bounds
figure('Name', 'Slip Rate Match', 'Position', [100 100 1100 500]);

% 3. Plot 95% Confidence Interval (p5-p95 shaded area)
fill([x_idx, fliplr(x_idx)], [p95_SR', fliplr(p5_SR')], ...
     [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
     'DisplayName', 'Ensemble 95% Bounds (p5-p95)'); 
hold on;

% 4. Plot Geological Targets (Max, Min, and Mean Target)
plot(targets.maxSR, '--', 'Color', [0 0.447 0.741], 'LineWidth', 1.2, 'DisplayName', 'Max Target');
plot(targets.minSR, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.2, 'DisplayName', 'Min Target');
meanTarget = (targets.maxSR + targets.minSR) / 2;
plot(meanTarget, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Mean Target');

% 5. Plot Ensemble Mean Model
plot(meanSR, 'r', 'LineWidth', 1.5, 'DisplayName', 'Model (Ensemble Mean)');

% Formatting
grid on; box on;
title('Slip Rate Inversion Results & Ensemble Stability', 'FontSize', 12);
ylabel('Slip Rate (m/yr)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Subsection Index', 'FontSize', 11, 'FontWeight', 'bold');
xlim([1, numSub]);
legend('Location', 'northeastoutside');

% --- Normalized Slip Rate Position Index (R_i Test) ---

% 1. Recupera i bounds geologici [numSub x 1]
SR_min = targets.minSR(:);
SR_max = targets.maxSR(:);
range_SR = SR_max - SR_min;

% Evita divisioni per zero nel caso improbabile in cui min == max
range_SR(range_SR == 0) = eps; 

% 2. Calcola la matrice R_i per TUTTE le 100 run [numSub x 100]
% all_runs_SR ha dimensione [numSub x 100]
R_all = (all_runs_SR - SR_min) ./ range_SR;

% 3. Calcola la media R_i per ogni sotto-sezione attraverso le 100 run [numSub x 1]
R_mean = mean(R_all, 2);

% --- VISUALIZZAZIONE ---
figure('Name', 'Normalized Slip Rate Position (R_i Analysis)', 'Position', [100 100 1100 450]);

% Pannello A: Istogramma della distribuzione di TUTTI i valori R_i (60 sub x 100 run)
subplot(1, 2, 1);
histogram(R_all(:), 20, 'FaceColor', [0.3 0.7 0.4], 'EdgeColor', 'w', 'Normalization', 'pdf');
hold on;
xline(0.5, 'r--', 'Mid-range (0.5)', 'LineWidth', 1.5, 'FontSize', 10, 'LabelVerticalAlignment', 'top');
grid on; box on;
title('a) Global Distribution of R_i Across Ensemble', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Normalized Slip Rate Index (R_i)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 10, 'FontWeight', 'bold');
xlim([-0.1 1.1]);

% Pannello B: Posizione media \overline{R_i} lungo il profilo delle sotto-sezioni
subplot(1, 2, 2);
plot(1:length(R_mean), R_mean, 'o-', 'Color', [0.1 0.4 0.6], 'LineWidth', 1.2, ...
     'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerSize', 5, 'DisplayName', '\overline{R_i} per Subsection');
hold on;
yline(0.5, 'r--', 'Ideal Center (0.5)', 'LineWidth', 1.5);
yline(0, 'k:', 'Min Bound (0)');
yline(1, 'k:', 'Max Bound (1)');
grid on; box on;
title('b) Mean Relative Position \overline{R_i} per Subsection', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Subsection Index', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Mean Relative Position \overline{R_i}', 'FontSize', 10, 'FontWeight', 'bold');
ylim([-0.1 1.1]);
xlim([1 length(R_mean)]);

sgtitle('Kinematic Allocation Unbiasedness Test Across Ensemble (N=100)', 'FontSize', 12, 'FontWeight', 'bold');
end
