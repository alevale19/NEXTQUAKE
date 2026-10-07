function plotMFDComparison(bestRates, invData, targets, params)
% PLOTMFDCOMPARISON - Plots Model MFD vs Target MFD with 5th and 95th Percentile Bounds.
% Correctly computes cumulative percentiles directly across ensemble realizations.

    fprintf('📊 Plotting MFD Comparison...\n');
    numBins = numel(targets.mfdBins);
    
    % --- 1. HANDLE ENSEMBLE INPUT (FULL MATRIX VS SINGLE VECTOR FALLBACK) ---
    if isstruct(bestRates) && isfield(bestRates, 'matrix')
        % ENSEMBLE MATRIX CASE: [numRuptures x numRealizations]
        ratesMatrix = bestRates.matrix;
        numRealizations = size(ratesMatrix, 2);
        
        % A. Compute incremental rates for each ensemble realization
        incMatrix = zeros(numBins, numRealizations);
        for k = 1:numBins
            binMask = (targets.rupBinIdx == k);
            incMatrix(k, :) = sum(ratesMatrix(binMask, :), 1);
        end
        
        % B. Compute cumulative rates for each individual realization
        cumMatrix = flipud(cumsum(flipud(incMatrix), 1));
        
        % C. Extract percentiles directly from the respective distributions
        modelInc_mean = mean(incMatrix, 2);
        modelInc_p5   = prctile(incMatrix, 5, 2);
        modelInc_p95  = prctile(incMatrix, 95, 2);
        
        modelCum_mean = mean(cumMatrix, 2);
        modelCum_p5   = prctile(cumMatrix, 5, 2);   % Native cumulative 5th percentile
        modelCum_p95  = prctile(cumMatrix, 95, 2);  % Native cumulative 95th percentile
        
        hasPercentiles = true;
        
    else
        % FALLBACK CASE: Single vector or simple mean structure without full ensemble
        if isstruct(bestRates)
            rates_mean = bestRates.mean;
        else
            rates_mean = bestRates;
        end
        
        modelInc_mean = zeros(numBins, 1);
        for k = 1:numBins
            binMask = (targets.rupBinIdx == k);
            modelInc_mean(k) = sum(rates_mean(binMask));
        end
        
        modelCum_mean = flipud(cumsum(flipud(modelInc_mean)));
        modelCum_p5   = modelCum_mean;
        modelCum_p95  = modelCum_mean;
        modelInc_p5   = modelInc_mean;
        modelInc_p95  = modelInc_mean;
        hasPercentiles = false;
    end
    
    % --- 2. CALCULATE TARGET CUMULATIVE MFD ---
    targetCumulative = flipud(cumsum(flipud(targets.mfdTarget)));
    
    % --- 3. GENERATE PLOTS ---
    figure('Name', 'MFD Comparison: Model vs Constraint', 'Color', 'w');
    
    % -----------------------------------------
    % SUBPLOT 1: INCREMENTAL RATES
    % -----------------------------------------
    subplot(2,1,1);
    bar(targets.mfdBins, modelInc_mean, 1, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', [0.4 0.4 0.4], 'DisplayName', 'Inversion (Mean)');
    hold on;
    
    if hasPercentiles
        errorbar(targets.mfdBins, modelInc_mean, modelInc_mean - modelInc_p5, modelInc_p95 - modelInc_mean, ...
            'k.', 'LineWidth', 1.2, 'CapSize', 4, 'DisplayName', '95% bounds');
    end
    
    plot(targets.mfdBins, targets.mfdTarget, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Regional MFD');
    set(gca, 'YScale', 'log'); grid on; box on;
    ylabel('Incremental Annual Rate [1/yr]');
    title('Incremental Magnitude-Frequency Distribution (MFD)');
    legend('Location', 'northeast');
    xlim([min(targets.mfdBins)-0.15, max(targets.mfdBins)+0.15]);
    ylim([1e-5, max([modelInc_p95; targets.mfdTarget])*5]);
    
    % -----------------------------------------
    % SUBPLOT 2: CUMULATIVE RATES (Corrected Bounds)
    % -----------------------------------------
    subplot(2,1,2);
    
    if hasPercentiles
        % Prepare polygon coordinates for the shaded uncertainty band
        xFill = [targets.mfdBins(:); flipud(targets.mfdBins(:))]';
        yFill = [modelCum_p5(:); flipud(modelCum_p95(:))]';
        
        % Plot percentile boundary lines
        plot(targets.mfdBins, modelCum_p5, ':', 'Color', [0.0 0.45 0.74], 'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(targets.mfdBins, modelCum_p95, ':', 'Color', [0.0 0.45 0.74], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    
    plot(targets.mfdBins, modelCum_mean, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Inversion (Mean)');
    hold on;
    fill(xFill, yFill, [0.0 0.45 0.74], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '95% bounds');
    plot(targets.mfdBins, targetCumulative, 'r--', 'LineWidth', 2, 'DisplayName', 'Regional MFD');
    
    set(gca, 'YScale', 'log'); grid on; box on;
    xlabel('Mw');
    ylabel('Cumulative Annual Rate [1/yr]');
    title('Cumulative Magnitude-Frequency Distribution (MFD)');
    legend('Location', 'northeast');
    xlim([min(targets.mfdBins)-0.15, max(targets.mfdBins)+0.15]);
    ylim([1e-5, max([modelCum_p95; targetCumulative])*5]);
    
    fprintf('✅ MFD Visualization complete.\n');
end
