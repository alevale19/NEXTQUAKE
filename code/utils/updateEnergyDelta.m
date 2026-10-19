function [deltaEnergy, newState] = updateEnergyDelta(r, deltaRate, oldRate, rates, state, invData, targets)
% UPDATEENERGYDELTA - Computes energy change in O(1) time updating only modified rows.
%   r         : Index of the modified rupture
%   deltaRate : Change in rate (newRate - oldRate)
%   oldRate   : Previous rate value rates(r)

    newState = state; % Copy current state
    newRate = oldRate + deltaRate;
    
    %% 1. DELTA SLIP RATE (Only for subsections involved in rupture r)
    % Find subsections intersected by rupture r
    subIndices = find(invData.Gsr(r, :) > 0);
    
    % Get old contribution to energy for these subsections
    old_sr = state.currentSR(subIndices);
    old_diff = zeros(size(old_sr));
    
    idx_b = old_sr < targets.minSR(subIndices);
    old_diff(idx_b) = (targets.minSR(subIndices(idx_b)) - old_sr(idx_b)) ./ targets.stdSR(subIndices(idx_b));
    idx_a = old_sr > targets.maxSR(subIndices);
    old_diff(idx_a) = (old_sr(idx_a) - targets.maxSR(subIndices(idx_a))) ./ targets.stdSR(subIndices(idx_a));
    
    % Update currentSR vector for affected subsections
    % deltaSR = Gsr(r, sub)^T * (deltaRate * Dr(r))
    addedSR = (invData.Gsr(r, subIndices)' .* invData.Dr(r)) * deltaRate;
    new_sr = old_sr + addedSR;
    
    % Calculate new contribution to energy for these subsections
    new_diff = zeros(size(new_sr));
    idx_b_new = new_sr < targets.minSR(subIndices);
    new_diff(idx_b_new) = (targets.minSR(subIndices(idx_b_new)) - new_sr(idx_b_new)) ./ targets.stdSR(subIndices(idx_b_new));
    idx_a_new = new_sr > targets.maxSR(subIndices);
    new_diff(idx_a_new) = (new_sr(idx_a_new) - targets.maxSR(subIndices(idx_a_new))) ./ targets.stdSR(subIndices(idx_a_new));
    
    % Delta Energy SR = sum(new_diff^2) - sum(old_diff^2)
    delta_energySR = sum(new_diff.^2) - sum(old_diff.^2);
    
    % Apply to state
    newState.currentSR(subIndices) = new_sr;
    newState.energySR = state.energySR + delta_energySR;

    %% 2. DELTA REGIONAL MFD (Only for the magnitude bin of rupture r)
    binIdx = targets.rupBinIdx(r);
    delta_energyMFD = 0;
    
    if ~isempty(targets.mfdTarget) && binIdx > 0
        old_bin_mfd = state.currentMFD(binIdx);
        new_bin_mfd = old_bin_mfd + deltaRate;
        
        old_mfd_diff = log10(old_bin_mfd + 1e-10) - log10(targets.mfdTarget(binIdx) + 1e-10);
        new_mfd_diff = log10(new_bin_mfd + 1e-10) - log10(targets.mfdTarget(binIdx) + 1e-10);
        
        delta_energyMFD = (new_mfd_diff^2 - old_mfd_diff^2) * targets.weightMFD;
        
        newState.currentMFD(binIdx) = new_bin_mfd;
        newState.energyMFD = state.energyMFD + delta_energyMFD;
    end

    %% 3. DELTA PALEOSEISMIC (Only for paleo sites intersected by rupture r)
    delta_energyPaleo = 0;
    if ~isempty(invData.Gpaleo)
        paleoSites = find(invData.Gpaleo(r, :) > 0);
        if ~isempty(paleoSites)
            old_paleo = state.currentPaleo(paleoSites);
            added_paleo = invData.Gpaleo(r, paleoSites)' * deltaRate;
            new_paleo = old_paleo + added_paleo;
            
            old_p_diff = (old_paleo - targets.paleoRate(paleoSites)) ./ targets.paleoStd(paleoSites);
            new_p_diff = (new_paleo - targets.paleoRate(paleoSites)) ./ targets.paleoStd(paleoSites);
            
            delta_energyPaleo = (sum(new_p_diff.^2) - sum(old_p_diff.^2)) * targets.weightPaleo;
            
            newState.currentPaleo(paleoSites) = new_paleo;
            newState.energyPaleo = state.energyPaleo + delta_energyPaleo;
        end
    end

    %% 4. DELTA SEGMENTATION
    delta_energySeg = 0;
    if ~isempty(invData.badRupIdx) && ismember(r, invData.badRupIdx)
        delta_energySeg = deltaRate * targets.weightSeg;
        newState.energySeg = state.energySeg + delta_energySeg;
    end

    %% 5. DELTA L2 REGULARIZATION
    % d/dx (x^2) -> (x + dx)^2 - x^2 = 2*x*dx + dx^2
    delta_energyL2 = targets.weightL2 * (2 * oldRate * deltaRate + deltaRate^2);
    newState.energyL2 = state.energyL2 + delta_energyL2;

    %% TOTAL DELTA ENERGY
    deltaEnergy = delta_energySR + delta_energyMFD + delta_energyPaleo + delta_energySeg + delta_energyL2;
    newState.totalEnergy = state.totalEnergy + deltaEnergy;
end