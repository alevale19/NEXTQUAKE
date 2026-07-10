function plotLocalMFD(bestRates, invData, targets, GSD, params)
% PLOTLOCALMFD - Calculates, formats, and saves standalone MFD plots 
% for each individual fault inside the '/Output/FaultMFD/' folder.
% Each plot compares the model's b-value slope against a Truncated GR target
% calibrated on the specific fault's geological moment rate.
    fprintf('📊 Generating and saving individual Fault MFD plots (Truncated GR Target)...\n');
    
    % --- Step A: Define and create output directory ---
    outputDir = './Output/FaultMFD';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
        fprintf('📁 Created directory: %s\n', outputDir);
    end
    
    % Constants from params
    shear_modulus = params.mu;
    MoRateReduction = params.MoRateReduction;
    
    allFaultNames = unique(GSD.FaultName);
    numFaults = numel(allFaultNames);
    numBins = numel(targets.mfdBins);
    
    % Hanks & Kanamori constants
    c_hk = 1.5; f_hk = 9.05;
    M_bins = 10.^(c_hk .* targets.mfdBins + f_hk);
    
    for f = 1:numFaults
        faultName = allFaultNames{f};
        idx = invData.faultRuptureIndices{f};
        
        % Filter GSD table rows for the current fault to compute its specific MoRate
        faultRows = strcmp(GSD.FaultName, faultName);
        
        % =========================================================
        % 1. COMPUTE FAULT SPECIFIC GEOLOGICAL MOMENT RATE (MoRate)
        % =========================================================
        faultMoRate = 0;
        if any(faultRows)
            L_m = GSD.Length_km(faultRows) * 1000;
            W_m = GSD.Width_km(faultRows) * 1000;
            
            % Check coupling property
            if ismember('Coupling', GSD.Properties.VariableNames)
                coupling = GSD.Coupling(faultRows);
                aseismicFactor = zeros(size(coupling));
                aseismicFactor(coupling < 1.0) = 0.4;
            else
                coupling = ones(size(L_m));
                aseismicFactor = zeros(size(L_m));
            end
            
            % Nominal slip rate from targets or GSD (ensure correct indexing)
            secSlipRate_m_yr = (targets.minSR(faultRows) + targets.maxSR(faultRows)) / 2;
            
            % Fault seismic moment rate
            secMoRate = shear_modulus .* (L_m .* (W_m .* (1 - aseismicFactor))) .* (secSlipRate_m_yr .* coupling);
            faultMoRate = sum(secMoRate) * (1 - MoRateReduction);
        end
        
        % =========================================================
        % 2. COMPUTE MODEL MFD FOR THIS FAULT
        % =========================================================
        faultIncremental = zeros(numBins, 1);
        minMag = Inf; maxMag = -Inf;
        
        if ~isempty(idx)
            f_mags = invData.Mw(idx);
            f_rates = bestRates.mean(idx);
            
            minMag = min(f_mags);
            maxMag = max(f_mags);
            
            for r = 1:numel(idx)
                [~, binIdx] = min(abs(targets.mfdBins - f_mags(r)));
                faultIncremental(binIdx) = faultIncremental(binIdx) + f_rates(r);
            end
        end
        faultCumulative = flipud(cumsum(flipud(faultIncremental)));
        
        % =========================================================
        % 3. GET TARGET B-VALUE & ESTIMATE CURRENT B-VALUE
        % =========================================================
        targetB = 1.0; % Default fallback
        currentB = NaN;
        
        if ~isempty(idx)
            subIdxs = invData.Gsr(idx(1),:) > 0;
            if isfield(targets, 'localBValue') && any(subIdxs)
                tB = targets.localBValue(subIdxs);
                if ~isempty(tB), targetB = tB(1); end
            end
            
            % Safety check for Aki MLE b-value calculation
            if sum(f_rates > 1e-8) > 2 && sum(f_rates) > 0
                currentB = estimateBValue(f_mags, f_rates);
            end
        end
        
        % =========================================================
        % 4. RECONSTRUCT THE TRUNCATED GR TARGET FOR THIS FAULT
        % =========================================================
        faultTargetCumulative = zeros(numBins, 1);
        
        if faultMoRate > 0 && ~isinf(minMag) && ~isinf(maxMag)
            % Compute scaling parameters bound to the specific fault's magnitude range
            Mt = 10^(c_hk * minMag + f_hk);
            mag_delta = targets.mfdBins(2) - targets.mfdBins(1);
            Mxp = 10^(c_hk * (maxMag + mag_delta) + f_hk); 
            
            Beta = (2/3) * targetB;
            
            % Core Truncated Gutenberg-Richter Formulation
            TruncGR = (((Mt ./ M_bins).^Beta) - ((Mt ./ Mxp).^Beta)) ./ (1 - (Mt ./ Mxp).^Beta);
            % Clean up values outside of valid range
            TruncGR(targets.mfdBins < minMag) = NaN;
            TruncGR(targets.mfdBins > maxMag) = 0;
            
            % Derive incremental steps to anchor total moment balancing
            TruncGR_valid = TruncGR; TruncGR_valid(isnan(TruncGR)) = 0;
            Incremental_GR = -diff([TruncGR_valid; 0]);
            
            Incremental_Morate = Incremental_GR .* M_bins;
            
            if sum(Incremental_Morate) > 0
                Incremental_Morate_balanced = (Incremental_Morate .* faultMoRate) ./ sum(Incremental_Morate);
                % Balanced target incremental rates
                faultTargetIncremental = (Incremental_GR .* Incremental_Morate_balanced) ./ (Incremental_Morate);
                faultTargetIncremental(isnan(faultTargetIncremental)) = 0;
                
                % Integrate to obtain cumulative target line
                faultTargetCumulative = flipud(cumsum(flipud(faultTargetIncremental)));
                faultTargetCumulative(targets.mfdBins < minMag | targets.mfdBins > maxMag) = NaN;
            end
        end
        
        % =========================================================
        % 5. CREATE STANDALONE FIGURE FOR THE FAULT
        % =========================================================
        hFig = figure('Name', faultName, 'Color', 'w', ...
                      'Position', [100 100 700 550], 'Visible', 'off');
        hold on;
        
        % Plot Truncated GR Target (Red dashed curve)
        h_targ = [];
        if any(faultTargetCumulative > 0, 'all')
            h_targ = plot(targets.mfdBins, faultTargetCumulative, 'r--', 'LineWidth', 1.8, ...
                'DisplayName', sprintf('Truncated GR Target (b = %.2f)', targetB));
        end
        
        % Plot Model Result (Blue solid line with filled circles)
        h_mod = plot(targets.mfdBins, faultCumulative, 'b-o', 'LineWidth', 2.0, ...
            'MarkerSize', 5, 'MarkerFaceColor', 'b', 'DisplayName', 'Model MFD');
        
        % Axis formatting
        set(gca, 'YScale', 'log', 'FontSize', 11);
        grid on; box on;
        
        maxVal = max([safeMax(faultCumulative), safeMax(faultTargetCumulative), 1e-3]);
        ylim([1e-5, maxVal * 5]);
        xlim([min(targets.mfdBins)-0.1, 8.0]);
      
        % Labels and Titles
        ylabel('Cumulative Rate (N \geq Mw)', 'FontSize', 12);
        xlabel('Magnitude (Mw)', 'FontSize', 12);
        
        if ~isnan(currentB)
            titleStr = {sprintf('Fault MFD: %s', faultName), ...
                        sprintf('b_{target} = %.2f, b_{model} = %.2f', targetB, currentB)};
        else
            titleStr = {sprintf('Fault MFD: %s', faultName), ...
                        sprintf('b_{target} = %.2f, b_{model} = NaN', targetB)};
        end
        
        title(titleStr, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'bold');
        
        % Build Legend
        if ~isempty(h_targ)
            legend([h_mod, h_targ], 'Location', 'southwest', 'FontSize', 11);
        else
            legend(h_mod, 'Location', 'southwest', 'FontSize', 11);
        end
        
        % =========================================================
        % 6. EXPORT AND CLEAN UP
        % =========================================================
        safeFaultName = regexprep(faultName, '[ \/\\:\*\?"<>\|]', '_');
        fileName = fullfile(outputDir, sprintf('MFD_%s.png', safeFaultName));
        
        saveas(hFig, fileName);
        close(hFig);
    end
    
    fprintf('✅ All individual fault MFD plots saved into "%s".\n', outputDir);
end

function b = estimateBValue(mags, rates)
    % Estimate b-value using weighted mean magnitude (Aki, 1965)
    mMin = min(mags);
    meanM = sum(rates .* mags) / sum(rates);
    b = (1 / (meanM - mMin)) * log10(exp(1));
    b = max(0.1, min(b, 2.0));
end

function m = safeMax(v)
    % Internal helper to handle arrays containing NaNs safely
    m = max(v(~isnan(v)));
    if isempty(m), m = 0; end
end