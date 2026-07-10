function [energy, components] = calculateEnergy(rates, invData, targets, params)
% CALCULATEENERGY - Advanced cost function with Box-Constraints, Paleo, 
% Local MFD (b-values), and Segmentation penalties.

    %% 1. SLIP RATE MISFIT (Box-Constraint Logic)
    % currentSR = [Subs x 1]
    currentSR = (invData.Gsr' * (rates .* invData.Dr));
    sr_diff = zeros(size(currentSR));
    
    % Penalty only if outside [minSR, maxSR]
    idx_below = currentSR < targets.minSR;
    sr_diff(idx_below) = (targets.minSR(idx_below) - currentSR(idx_below)) ./ targets.stdSR(idx_below);
    
    idx_above = currentSR > targets.maxSR;
    sr_diff(idx_above) = (currentSR(idx_above) - targets.maxSR(idx_above)) ./ targets.stdSR(idx_above);
    
    energySR = sum(sr_diff.^2);

    %% 2. REGIONAL MFD MISFIT (Log-Linear)
    currentMFD = zeros(size(targets.mfdBins));
    for i = 1:numel(targets.mfdBins)
        currentMFD(i) = sum(rates(targets.rupBinIdx == i));
    end
    
    if ~isempty(targets.mfdTarget)
        mfd_diff = (log10(currentMFD + 1e-10) - log10(targets.mfdTarget + 1e-10));
        energyMFD = sum(mfd_diff.^2) * targets.weightMFD;
    else
        energyMFD = 0;
    end

    %% 3. PALEOSEISMIC MISFIT
    % Compare predicted event rates with observed paleo-rates at specific sites
    if ~isempty(invData.Gpaleo)
        % currentPaleo = [NumSites x 1]
        currentPaleo = invData.Gpaleo' * rates; 
        paleo_diff = (currentPaleo - targets.paleoRate) ./ targets.paleoStd;
        energyPaleo = sum(paleo_diff.^2) * targets.weightPaleo;
    else
        energyPaleo = 0;
    end

%% 4. LOCAL MFD MISFIT (b-value control per fault)
    % Penalize the deviation from the target b-value for each fault section
    energyLocalMFD = 0;
    if targets.weightLocalMFD > 0
        for f = 1:numel(invData.faultRuptureIndices)
            idx = invData.faultRuptureIndices{f};
            if isempty(idx), continue; end

            f_rates = rates(idx);
            f_mags = invData.Mw(idx);

            % Target b-value for this specific fault (from first subsection involved)
            % We use the first index as b-values are constant per fault
            targetB = targets.localBValue(invData.Gsr(idx(1),:) > 0); 
            if isempty(targetB), targetB = 1.0; else, targetB = targetB(1); end

            % Calculate current local b-value (simplified MLE or slope)
            % If fewer than 2 ruptures are active, we skip to avoid noise
            if sum(f_rates > 1e-8) > 2
                currentB = estimateBValue(f_mags, f_rates);
                energyLocalMFD = energyLocalMFD + (currentB - targetB)^2;
            end
        end
        energyLocalMFD = energyLocalMFD * targets.weightLocalMFD;
    end

    %% 5. SEGMENTATION PENALTY
    % Directly penalize the rates of "forbidden" ruptures
    if ~isempty(invData.badRupIdx)
        % Linear penalty: the more slip on bad ruptures, the higher the energy
        energySeg = sum(rates(invData.badRupIdx)) * targets.weightSeg;
    else
        energySeg = 0;
    end

    %% 6. L2 REGULARIZATION (SMOOTHING)
    energyL2 = targets.weightL2 * sum(rates.^2);

    %% 7. TOTAL ENERGY
    energy = energySR + energyMFD + energyPaleo + energyLocalMFD + energySeg + energyL2;

    %% 8. DIAGNOSTICS
    components.SR = energySR;
    components.MFD = energyMFD;
    components.Paleo = energyPaleo;
    components.LocalMFD = energyLocalMFD;
    components.Seg = energySeg;
    components.L2 = energyL2;
end

function b = estimateBValue(mags, rates)
    % Estimate b-value using weighted mean magnitude (Aki, 1965)
    mMin = min(mags);
    meanM = sum(rates .* mags) / sum(rates);
    b = (1 / (meanM - mMin)) * log10(exp(1));
    % Cap b-value to avoid numerical explosions
    b = max(0.1, min(b, 2.0));
end