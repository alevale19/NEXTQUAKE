function [bestRates, energyHistory] = runSimulatedAnnealing(invData, targets, params, faultResults)
% RUNSIMULATEDANNEALING - Parallel Ensemble with Fast Delta-Energy Updating
% Minimizes total energy using O(1) row-based energy updates for high-performance SA.

    % --- 1. HYPERPARAMETERS & SEED SETUP ---
    numIter = params.numIter;
    numRuns = params.runSA;
    T_init  = 10.0;
    T_final = 0.001;
    cooling_rate = (T_final / T_init)^(1 / numIter); 
    
    numRup = invData.numRuptures;
    numFaults = length(faultResults);
    
    % --- SEISMOTECTONIC STATIC DELTASCALE CALIBRATION ---
    validSR = targets.maxSR(targets.maxSR > 0);
    if ~isempty(validSR)
        referenceSR = prctile(validSR, 10); 
        staticDeltaScale = referenceSR; 
    else
        staticDeltaScale = 1e-4; 
    end
    
    % Filtering valid ruptures (exclude segmentation violations)
    validRupIdx = setdiff(1:numRup, invData.badRupIdx);
    numValid = numel(validRupIdx);
    
    % =========================================================================
    % ⚠️ INTELLIGENT WARNING: INITIAL ITERATION SUFFICIENCY CHECK
    % =========================================================================
    samplingFactorK = numIter / numValid;
    minRecommendedK = params.minRecommendedK;
    
    if samplingFactorK < minRecommendedK
        warning('⚠️ INSUFFICIENT ITERATIONS DETECTED relative to the scale of the problem!');
        fprintf('   ---------------------------------------------------------------------------\n');
        fprintf('   👉 You have set %d iterations for %d valid ruptures.\n', numIter, numValid);
        fprintf('   👉 Each rupture will be sampled an average of %.1f times.\n', samplingFactorK);
        fprintf('   🔔 SUGGESTION: Increase the number of iterations to at least %d for this configuration,\n', numValid * 2000);
        fprintf('      or proceed knowing that the results might be underfitted.\n');
        fprintf('   ---------------------------------------------------------------------------\n\n');
        pause(3); 
    end

    % Pre-allocate ensemble variables for parfor
    ensembleRates = zeros(numRup, numRuns); 
    energyHistory = zeros(numIter, numRuns); 
    faultSlipEnsemble = zeros(numFaults, numRuns);
    
    fprintf('\n🚀 Starting Parallel SA Ensemble (%d Runs)...\n', numRuns);
    
   % --- 2. PARALLEL ENSEMBLE LOOP ---
    totalAnnealingTimer = tic; % Timer globale dell'ensemble
    
    parfor n = 1:numRuns
        fprintf('▶️ Worker: starting Run %d...\n', n);
        
        if isfield(params, 'seed')
            rng(params.seed + n);
        end
        
        localRates = zeros(numRup, 1);
        state = initEnergyState(localRates, invData, targets);
        
        bestRatesInRun = localRates;
        bestEnergyInRun = state.totalEnergy;
        T = T_init;
        
        % --- TIMER WORKER LOCALE ---
        runTimer = tic; 
        
        % --- Main Optimization Loop ---
        for i = 1:numIter
            idx = validRupIdx(randi(numValid));
            oldRate = localRates(idx);
            delta = (rand() - 0.5) * staticDeltaScale;
            
            if oldRate + delta < 0
                delta = -oldRate;
            end
            
            [deltaE, candidateState] = updateEnergyDelta(idx, delta, oldRate, localRates, state, invData, targets);
            
            if deltaE < 0 || rand() < exp(-deltaE / T)
                localRates(idx) = oldRate + delta;
                state = candidateState;
                
                if state.totalEnergy < bestEnergyInRun
                    bestEnergyInRun = state.totalEnergy;
                    bestRatesInRun = localRates;
                end
            end
            
            T = T * cooling_rate;
            energyHistory(i, n) = state.totalEnergy; 
        end
        
        % --- CALCOLO VELOCITÀ WORKER ---
        runTime = toc(runTimer);
        workerItsPerSec = numIter / runTime;
        
        ensembleRates(:, n) = bestRatesInRun; 
        
        % Final Fault-Level Slip Rate Calculation
        currentSR_sub = (invData.Gsr' * (bestRatesInRun .* invData.Dr));
        localFaultSlip = zeros(numFaults, 1);
        for fIdx = 1:numFaults
            subIdx = faultResults(fIdx).GlobalSubIndices;
            if ~isempty(subIdx)
                localFaultSlip(fIdx) = mean(currentSR_sub(subIdx));
            end
        end
        faultSlipEnsemble(:, n) = localFaultSlip; 
        
        fprintf('✅ Worker: Run %d finished in %.2f s (⚡ %.0f iter/sec per core). Final Energy: %.2e\n', ...
                n, runTime, workerItsPerSec, bestEnergyInRun);
    end
    
    % --- REPORT PRESTAZIONI GLOBALI CLUSTER ---
    totalTime = toc(totalAnnealingTimer);
    totalIterationsAllRuns = numIter * numRuns;
    globalItsPerSec = totalIterationsAllRuns / totalTime;
    
    fprintf('\n=========================================================\n');
    fprintf('⚡ PERFORMANCE SUMMARY (Simulated Annealing Engine):\n');
    fprintf('   Total Runs Executed:       %d\n', numRuns);
    fprintf('   Total Iterations Processed: %d (%d per run)\n', totalIterationsAllRuns, numIter);
    fprintf('   Total Execution Time:      %.2f seconds\n', totalTime);
    fprintf('   🚀 Aggregated Throughput:   %.0f iter/sec (across all parallel workers)\n', globalItsPerSec);
    fprintf('=========================================================\n\n');
    
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
    
    fprintf('🏆 Global Minimum Solution Found in Run #%d.\n', absoluteBestRunIdx);

    % --- 4. DIAGNOSTIC PLOTS & STABILITY ANALYSIS ---
    mean_rup_all = mean(ensembleRates, 2);
    std_rup_all  = std(ensembleRates, 0, 2);
    active_rup   = mean_rup_all > 1e-6; 
    cv_rup_vec   = std_rup_all(active_rup) ./ mean_rup_all(active_rup); 

    all_runs_SR  = invData.Gsr' * (ensembleRates .* invData.Dr);
    mean_fault   = mean(all_runs_SR, 2);
    std_fault    = std(all_runs_SR, 0, 2);
    active_fault = mean_fault > 1e-6;
    cv_fault_vec = std_fault(active_fault) ./ mean_fault(active_fault); 

    % Figure 1: Ensemble CV Histograms
    figure('Name', 'Ensemble CV Histograms', 'Color', 'w');
    subplot(1, 2, 1);
    histogram(cv_rup_vec, 30, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w');
    grid on; box on;
    title('a) Rupture-Level Rate Variability', 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Coefficient of Variation (CV)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Number of Ruptures', 'FontSize', 10, 'FontWeight', 'bold');
    if ~isempty(cv_rup_vec)
        xline(median(cv_rup_vec), 'r--', sprintf(' Median: %.2f', median(cv_rup_vec)), ...
              'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', ...
              'LabelHorizontalAlignment', 'left', 'FontSize', 10);
    end

    subplot(1, 2, 2);
    histogram(cv_fault_vec, 15, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'w');
    grid on; box on;
    title('b) Subsection-Level Slip Rate Variability', 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Coefficient of Variation (CV)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Number of Subsections', 'FontSize', 10, 'FontWeight', 'bold');
    if ~isempty(cv_fault_vec)
        xline(median(cv_fault_vec), 'r--', sprintf(' Median: %.2f', median(cv_fault_vec)), ...
              'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', ...
              'LabelHorizontalAlignment', 'left', 'FontSize', 10);
    end
    sgtitle('Ensemble Solution Variability Across Independent Runs', 'FontSize', 12, 'FontWeight', 'bold');

    % Figure 2: Slip Rate Match with Bounds
    meanSR = mean(all_runs_SR, 2);
    p5_SR  = prctile(all_runs_SR, 5, 2);
    p95_SR = prctile(all_runs_SR, 95, 2);
    numSub = length(meanSR);
    x_idx  = 1:numSub;

    figure('Name', 'Slip Rate Match', 'Color', 'w');
    fill([x_idx, fliplr(x_idx)], [p95_SR', fliplr(p5_SR')], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
         'DisplayName', 'Ensemble 95% Bounds (p5-p95)'); 
    hold on;
    plot(targets.maxSR, '--', 'Color', [0 0.447 0.741], 'LineWidth', 1.2, 'DisplayName', 'Max Target');
    plot(targets.minSR, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.2, 'DisplayName', 'Min Target');
    meanTarget = (targets.maxSR + targets.minSR) / 2;
    plot(meanTarget, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Mean Target');
    plot(meanSR, 'r', 'LineWidth', 1.5, 'DisplayName', 'Model (Ensemble Mean)');
    grid on; box on;
    title('Slip Rate Inversion Results & Ensemble Stability', 'FontSize', 12);
    ylabel('Slip Rate (m/yr)', 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Subsection Index', 'FontSize', 11, 'FontWeight', 'bold');
    xlim([1, numSub]);
    legend('Location', 'northeastoutside');

    % Figure 3: Normalized Slip Rate Position Index (R_i Test)
    SR_min = targets.minSR(:);
    SR_max = targets.maxSR(:);
    range_SR = SR_max - SR_min;
    range_SR(range_SR == 0) = eps; 

    R_all = (all_runs_SR - SR_min) ./ range_SR;
    R_mean = mean(R_all, 2);

    figure('Name', 'Normalized Slip Rate Position (Ri Analysis)', 'Color', 'w');
    plot(1:length(R_mean), R_mean, 'o-', 'Color', [0.1 0.4 0.6], 'LineWidth', 1.2, ...
         'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerSize', 5, 'DisplayName', '\overline{R_i} per Subsection');
    hold on;
    yline(0, 'k:', 'Min Bound (0)');
    yline(1, 'k:', 'Max Bound (1)');
    grid on; box on;
    title('Mean Relative Position of R_i per Subsection', 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Subsection Index', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('R_i', 'FontSize', 10, 'FontWeight', 'bold');
    ylim([-0.1 1.1]);
    xlim([1 length(R_mean)]);
end
