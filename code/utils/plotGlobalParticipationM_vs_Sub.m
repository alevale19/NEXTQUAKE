function plotGlobalParticipationM_vs_Sub(outputs, invData, bestRates, subsectionData, isTimeDependent, params, targets)
% PLOTGLOBALPARTICIPATIONM_VS_SUB - Generates UCERF3-style global system profiles
% supporting both Poisson (Stationary) and BPT (Time-Dependent) probability modulations.
%
% Inputs:
%   outputs        : Structure from processAnnealingResults containing Subsections table
%   invData        : Structure containing rupture metadata (Mw, Gsr)
%   bestRates      : Vector of optimized rupture rates from Simulated Annealing
%   subsectionData : Table containing global subsections data
%   isTimeDependent: Logical flag. True for BPT-scaled rates, False for Poisson rates

    if nargin < 5, isTimeDependent = false; end

    if isTimeDependent
        fprintf('📊 Generating Global Magnitude vs Subsection Participation Plots (BPT TIME-DEPENDENT)...\n');
        plotTitleStr = 'BPT Time-Dependent';
    else
        fprintf('📊 Generating Global Magnitude vs Subsection Participation Plots (POISSON STATIONARY)...\n');
        plotTitleStr = 'Poisson Stationary';
    end

    numSubsections = height(subsectionData);
    numRuptures = length(bestRates.mean);
    logicalGsr = invData.Gsr > 0; 
    activeRupIdx = find(bestRates.mean > 0);
    
    % --- Define Magnitude bins (from 6.0 to 8.2 with 0.1 step) ---
    minmag = min(targets.mfdBins)-0.15;
    maxmag = max(targets.mfdBins)+0.15;
    magBins = minmag:params.mag_delta:maxmag;
    numMagBins = length(magBins);
    
    % --- Initialize the Incremental Rate Matrix (Magnitude x Subsections) ---
    incrementalMatrix = zeros(numMagBins, numSubsections);
    
    % --- Step 1: Populate the Incremental Matrix ---
    for r = activeRupIdx'
        rate = bestRates.mean(r);
        subsInvolved = find(logicalGsr(r, :));
        if isempty(subsInvolved), continue; end
        
        % SE TIME-DEPENDENT: Scaliamo il tasso originario con il Probability Gain
        if isTimeDependent
            P_poisson = outputs.Subsections.PoissonProb(subsInvolved);
            P_bpt     = outputs.Subsections.BPT_TimeDependentProb(subsInvolved);
            
            % Calcoliamo il Gain cella per cella gestendo in sicurezza le divisioni per zero
            G_subs = zeros(size(P_poisson));
            validGain = P_poisson > 0;
            G_subs(validGain) = P_bpt(validGain) ./ P_poisson(validGain);
            
            % Applichiamo il guadagno medio delle sottosezioni coinvolte al tasso della rottura
            rate = rate * mean(G_subs);
        end
        
        mw = invData.Mw(r);
        % Find which magnitude bin this rupture belongs to (nearest neighbor approach)
        [~, magIdx] = min(abs(magBins - mw));
        
        % Accumulate the modulated rate to all involved subsections
        incrementalMatrix(magIdx, subsInvolved) = incrementalMatrix(magIdx, subsInvolved) + rate;
    end
    
    % --- Step 2: Compute the Cumulative Matrix (Summing from high Mw down to low Mw) ---
    cumulativeMatrix = zeros(numMagBins, numSubsections);
    for m = 1:numMagBins
        cumulativeMatrix(m, :) = sum(incrementalMatrix(m:end, :), 1);
    end
    
    % --- Step 3: Convert to Log10 and manage zeros/extremely low values safely ---
    logIncremental = log10(incrementalMatrix);
    logIncremental(incrementalMatrix == 0 | logIncremental < -6) = -6; % Floor at 10^-6
    logCumulative = log10(cumulativeMatrix);
    logCumulative(cumulativeMatrix == 0 | logCumulative < -6) = -6;    % Floor at 10^-6
    
    % =====================================================================
    % --- RENDER THE UCERF3 MASTER PLOT ---
    % =====================================================================
    figName = sprintf('NextQuake - Global System Participation Profiles (%s)', plotTitleStr);
    figure('Name', figName, 'Color', 'w', 'Position', [100, 100, 900, 600]);
    
    % --- TOP PANEL: Incremental Rate ---
    ax1 = subplot(2, 1, 1);
    imagesc(1:numSubsections, magBins, logIncremental);
    set(gca, 'YDir', 'normal'); 
    grid on; 
    colormap(gca, jet); 
    clim([-6 -1]);      
    title(sprintf('Incremental rate (%s)', plotTitleStr), 'FontWeight', 'bold', 'FontSize', 11, 'HorizontalAlignment', 'center');
    ylabel('Magnitude');
    
    % Draw vertical dividers to separate parent faults visually
    hold on;
    faultTransitions = find(diff(subsectionData.Fault_ID) ~= 0);
    for t = 1:length(faultTransitions)
        plot([faultTransitions(t) faultTransitions(t)], [magBins(1) magBins(end)], 'w-', 'LineWidth', 0.8);
    end
    
    % --- BOTTOM PANEL: Cumulative Rate ---
    ax2 = subplot(2, 1, 2);
    imagesc(1:numSubsections, magBins, logCumulative);
    set(gca, 'YDir', 'normal');
    grid on;
    colormap(gca, jet);
    clim([-6 -1]);
    title(sprintf('Cumulative rate (%s)', plotTitleStr), 'FontWeight', 'bold', 'FontSize', 11, 'HorizontalAlignment', 'center');
    xlabel('Subsection ID'); ylabel('Magnitude');
    
    % Draw fault vertical dividers here too
    hold on;
    for t = 1:length(faultTransitions)
        plot([faultTransitions(t) faultTransitions(t)], [magBins(1) magBins(end)], 'w-', 'LineWidth', 0.8);
    end
    
    % --- Add the unified colorbar at the bottom ---
    cb = colorbar(ax2, 'eastoutside');
    ylabel(cb, 'Log_{10} participation rate (events/yr)', 'FontSize', 10, 'FontWeight', 'bold');

    ax1.Position(3) = ax2.Position(3);
    
    fprintf('✅ %s Global profile plotted successfully.\n', plotTitleStr);
end