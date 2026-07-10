function plotMFDComparison(bestRates, invData, targets, params)
% PLOTMFDCOMPARISON - Plots Model MFD vs Target MFD with 5th and 95th Percentile Bounds.
% Supports both vector input (classic) and structured ensemble input (percentiles).

    fprintf('📊 Plotting MFD Comparison with Epistemic Uncertainty...\n');
    numBins = numel(targets.mfdBins);
    
    % --- 1. HANDLE INPUT TYPE (STRUCTURED ENSEMBLE VS CLASSIC VECTOR) ---
    % If the user passes a simple vector, we fallback safely to using it as the mean
    if isstruct(bestRates)
        rates_mean = bestRates.mean;
        rates_p5   = bestRates.p5;
        rates_p95  = bestRates.p95;
        hasPercentiles = true;
    else
        rates_mean = bestRates;
        rates_p5   = bestRates; % Fallback
        rates_p95  = bestRates; % Fallback
        hasPercentiles = false;
    end

    % --- 2. CALCULATE INCREMENTAL RATES FOR EACH BOUND ---
    modelInc_mean = zeros(numBins, 1);
    modelInc_p5   = zeros(numBins, 1);
    modelInc_p95  = zeros(numBins, 1);
    
    for k = 1:numBins
        binMask = (targets.rupBinIdx == k);
        modelInc_mean(k) = sum(rates_mean(binMask));
        modelInc_p5(k)   = sum(rates_p5(binMask));
        modelInc_p95(k)  = sum(rates_p95(binMask));
    end
    
    % --- 3. CALCULATE CUMULATIVE RATES (N >= M) ---
    modelCum_mean = flipud(cumsum(flipud(modelInc_mean)));
    modelCum_p5   = flipud(cumsum(flipud(modelInc_p5)));
    modelCum_p95  = flipud(cumsum(flipud(modelInc_p95)));
    targetCumulative = flipud(cumsum(flipud(targets.mfdTarget)));
    
    % --- 4. GENERATE PLOTS ---
    figure('Name', 'MFD Comparison: Model Ensemble vs Target', 'Color', 'w', 'Position', [100 100 850 700]);
    
    % -----------------------------------------
    % SUBPLOT 1: INCREMENTAL RATES
    % -----------------------------------------
    subplot(2,1,1);
    
    % Plot the Mean Model as a baseline gray bar chart
    bar(targets.mfdBins, modelInc_mean, 1, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', [0.4 0.4 0.4], 'DisplayName', 'Model (Mean)');
    hold on;
    
    % If percentiles exist, plot them as uncertainty error whiskers or lines
    if hasPercentiles
        errorbar(targets.mfdBins, modelInc_mean, modelInc_mean - modelInc_p5, modelInc_p95 - modelInc_mean, ...
            'k.', 'LineWidth', 1.2, 'CapSize', 4, 'DisplayName', 'Model (5th-95th %ile)');
    end
    
    % Plot Target
    plot(targets.mfdBins, targets.mfdTarget, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Target (Geological/GR)');
    
    set(gca, 'YScale', 'log');
    grid on; box on;
    ylabel('Annual Rate (Incremental) [1/yr]');
    title('Incremental Magnitude-Frequency Distribution (MFD)');
    legend('Location', 'northeast');
    xlim([min(targets.mfdBins)-0.15, max(targets.mfdBins)+0.15]);
    ylim([1e-5, max([modelInc_p95; targets.mfdTarget])*5]);

    % -----------------------------------------
    % SUBPLOT 2: CUMULATIVE RATES (The classic Gutenberg-Richter Plot)
    % -----------------------------------------
    subplot(2,1,2);
    
    % If percentiles exist, draw a beautiful shaded uncertainty band
    if hasPercentiles
        % Prepare X and Y vectors for the fill polygon
        xFill = [targets.mfdBins(:); flipud(targets.mfdBins(:))]';
        yFill = [modelCum_p5(:); flipud(modelCum_p95(:))]';
        
        % Shaded blue area with transparency (Alpha = 0.15)
        fill(xFill, yFill, [0.0 0.45 0.74], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Model Epistemic Range (p5-p95)');
        hold on;
        
        % Plot individual percentile boundary lines
        plot(targets.mfdBins, modelCum_p5, ':', 'Color', [0.0 0.45 0.74], 'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(targets.mfdBins, modelCum_p95, ':', 'Color', [0.0 0.45 0.74], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    
    % Plot Mean Model and Target Cumulative Curves
    plot(targets.mfdBins, modelCum_mean, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Model (Ensemble Mean)');
    hold on;
    plot(targets.mfdBins, targetCumulative, 'r--', 'LineWidth', 2, 'DisplayName', 'Target (Cumulative)');
    
    set(gca, 'YScale', 'log');
    grid on; box on;
    xlabel('Magnitude (Mw)');
    ylabel('Annual Rate (N \geq Mw) [1/yr]');
    title('Cumulative Magnitude-Frequency Distribution (MFD)');
    legend('Location', 'northeast');
    xlim([min(targets.mfdBins)-0.15, max(targets.mfdBins)+0.15]);
    
    % Set a clean, dynamic Y-limit clamped to a minimum of 10^-5
    ylim([1e-5, max(targetCumulative)*5]);
    
    fprintf('✅ MFD Visualization with percentiles complete.\n');
end