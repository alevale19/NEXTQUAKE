function state = initEnergyState(rates, invData, targets)
% INITENERGYSTATE - Computes full initial vectors and energy components.

    % 1. Slip Rate
    state.currentSR = invData.Gsr' * (rates .* invData.Dr);
    sr_diff = zeros(size(state.currentSR));
    idx_below = state.currentSR < targets.minSR;
    sr_diff(idx_below) = (targets.minSR(idx_below) - state.currentSR(idx_below)) ./ targets.stdSR(idx_below);
    idx_above = state.currentSR > targets.maxSR;
    sr_diff(idx_above) = (state.currentSR(idx_above) - targets.maxSR(idx_above)) ./ targets.stdSR(idx_above);
    state.energySR = sum(sr_diff.^2);

    % 2. Regional MFD
    state.currentMFD = zeros(size(targets.mfdBins));
    for i = 1:numel(targets.mfdBins)
        state.currentMFD(i) = sum(rates(targets.rupBinIdx == i));
    end
    if ~isempty(targets.mfdTarget)
        mfd_diff = (log10(state.currentMFD + 1e-10) - log10(targets.mfdTarget + 1e-10));
        state.energyMFD = sum(mfd_diff.^2) * targets.weightMFD;
    else
        state.energyMFD = 0;
    end

    % 3. Paleoseismic
    if ~isempty(invData.Gpaleo)
        state.currentPaleo = invData.Gpaleo' * rates; 
        paleo_diff = (state.currentPaleo - targets.paleoRate) ./ targets.paleoStd;
        state.energyPaleo = sum(paleo_diff.^2) * targets.weightPaleo;
    else
        state.currentPaleo = [];
        state.energyPaleo = 0;
    end

    % 4. Segmentation
    if ~isempty(invData.badRupIdx)
        state.energySeg = sum(rates(invData.badRupIdx)) * targets.weightSeg;
    else
        state.energySeg = 0;
    end

    % 5. L2 Regularization
    state.energyL2 = targets.weightL2 * sum(rates.^2);

    % Total initial energy
    state.totalEnergy = state.energySR + state.energyMFD + state.energyPaleo + ...
                        state.energySeg + state.energyL2;
end