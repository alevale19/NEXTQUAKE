function plotInversionAnalysisMaps(outputs, faultResults, subsectionData, invData, bestRates, params)
% PLOTINVERSIONANALYSISMAPS - Advanced spatial hazard visualizations optimized for 3D meshes:
% 1. Model-Derived Fault Slip Rates (Parent Fault Level)
% 2. Mean Return Periods (Parent Fault Level)
% 3. Poisson Probabilities (Parent Fault Level)
% 4. BPT Time-Dependent Probabilities (Parent Fault Level)
% 5. Probability Gain Map (3D Subsection Patch Level)
% 6. Participation Rate Map (3D Subsection Patch Level in Log10 Scale)

    % --- Extract Parameters from Struct or Set Defaults ---
    if nargin < 6, params = struct(); end
    if isfield(params, 'forTwind'), forecastWindow = params.forTwind; else, forecastWindow = 30; end
    if isfield(params, 'forecastMw'), Mw_limit = params.forecastMw; else, Mw_limit = 6.5; end 
    
    mu = params.mu;
    numFaults = length(faultResults);
    numSubsections = height(subsectionData);
    
    % --- Extract Pre-Calculated Values from the Unified Output Engine ---
    faultParticipationRate = outputs.Faults.ParticipationRate;
    faultPoissonProb       = outputs.Faults.PoissonProb;
    faultBPTProb           = outputs.Faults.BPT_TimeDependentProb;
    
    subPoissonProb         = outputs.Subsections.PoissonProb;
    subBPTProb             = outputs.Subsections.BPT_TimeDependentProb;
    subParticipationRate   = outputs.Subsections.ParticipationRate;
    
    % --- Calculate Parent Fault Mean Return Periods ---
    returnPeriods = 1 ./ faultParticipationRate;
    returnPeriods(faultParticipationRate == 0) = NaN;
    
    % --- Rigorous Pre-Calculation of Subsection Probability Gain ---
    subGain = zeros(numSubsections, 1);
    validGainIdx = subPoissonProb > 0;
    subGain(validGainIdx) = subBPTProb(validGainIdx) ./ subPoissonProb(validGainIdx);
    
    % --- Back-Calculate Fault Slip Rates from Inversion Moment Rates ---
    faultModelSlipRate = zeros(numFaults, 1);
    activeRupIdx = find(bestRates.mean > 0);
    logicalGsr = invData.Gsr > 0; 
    for i = 1:numFaults
        currentFault = faultResults(i); 
        subIdx = currentFault.GlobalSubIndices;
        if isempty(subIdx), continue; end
        
        fArea = currentFault.TotalArea_km2; 
        totalFaultMoRate = 0;
        
        for r = activeRupIdx'
            if any(logicalGsr(r, subIdx))
                totalFaultMoRate = totalFaultMoRate + bestRates.mean(r) * invData.Dr(r) * fArea * mu;
            end
        end
        faultModelSlipRate(i) = totalFaultMoRate / (fArea * mu);
    end
    
    fprintf('📊 Generating Spatial Analysis Maps (Mw >= %.1f, Window = %d years)...\n', Mw_limit, forecastWindow);    
    
    % =====================================================================
    % MAP 1: MODEL SLIP RATE (Parent Fault Level)
    % =====================================================================
    figure('Name', 'NextQuake - Model Slip Rate', 'Color', 'w');
    hold on; grid on; axis equal; colormap(flipud(hot));
    for i = 1:numFaults
        f = faultResults(i);
        c = get_color_local(faultModelSlipRate(i), faultModelSlipRate, flipud(hot));
        plot(f.coordsWGS(:,1), f.coordsWGS(:,2), '-', 'LineWidth', 3, 'Color', c);
        text(f.coordsWGS(1,1), f.coordsWGS(1,2), f.Name, 'FontSize', 7, 'Interpreter', 'none');
    end
    title('Model Slip Rates (m/yr) - Inversion Result'); 
    cb = colorbar; ylabel(cb, 'Slip Rate [m/yr]');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
    % =====================================================================
    % MAP 2: RECURRENCE INTERVAL / RETURN PERIOD (Parent Fault Level)
    % =====================================================================
    figure('Name', 'NextQuake - Return Period', 'Color', 'w');
    hold on; grid on; axis equal; colormap(jet);
    for i = 1:numFaults
        f = faultResults(i);
        if ~isnan(returnPeriods(i)) && returnPeriods(i) < inf
            c = get_color_local(returnPeriods(i), returnPeriods(~isnan(returnPeriods) & returnPeriods < inf), jet);
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), '-', 'LineWidth', 3, 'Color', c);
            text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf('TR: %.0f y', returnPeriods(i)), 'FontSize', 7);
        else
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), ':', 'Color', [0.7 0.7 0.7]);
        end
    end
    title(sprintf('Mean Return Period (Years) for Mw \\geq %.1f', Mw_limit));
    cb = colorbar; set(gca, 'ColorScale', 'log'); ylabel(cb, 'Years');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
    % =====================================================================
    % MAP 3: POISSON PROBABILITY (Parent Fault Level)
    % =====================================================================
    maxP = max(faultPoissonProb(faultParticipationRate > 0));
    figure('Name', sprintf('NextQuake - Poisson Probability %dy', forecastWindow), 'Color', 'w');
    hold on; grid on; axis equal; colormap(parula);
    for i = 1:numFaults
        f = faultResults(i);
        if faultParticipationRate(i) > 0
            c = get_color_local(faultPoissonProb(i), faultPoissonProb(faultParticipationRate > 0), parula);
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), '-', 'LineWidth', 3.5, 'Color', c);
            text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf('P_p: %.1f%%', faultPoissonProb(i)*100), ...
                'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'none');
        else
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), ':', 'Color', [0.7 0.7 0.7]);
        end
    end
    title(sprintf('%d-Year Poisson Probability (Mw \\geq %.1f)', forecastWindow, Mw_limit));
    
    clim([0 maxP]);
    cb = colorbar; ylabel(cb, 'Probability');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
    % =====================================================================
    % MAP 4: BPT TIME-DEPENDENT PROBABILITY (Parent Fault Level)
    % =====================================================================
    maxPbpt = max(faultBPTProb(faultParticipationRate > 0));
    figure('Name', sprintf('NextQuake - BPT Time-Dependent Probability %dy', forecastWindow), 'Color', 'w');
    hold on; grid on; axis equal; colormap(spring);
    
    for i = 1:numFaults
        f = faultResults(i);
        if faultParticipationRate(i) > 0
            c = get_color_local(faultBPTProb(i), faultBPTProb(faultParticipationRate > 0), spring);
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), '-', 'LineWidth', 3.5, 'Color', c);
            text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf('%s (BPT: %.1f%%)', f.Name, faultBPTProb(i)*100), ...
                'FontSize', 7, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1], 'Interpreter', 'none');
        else
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), ':', 'Color', [0.7 0.7 0.7]);
            text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf('%s (BPT: 0.0%%)', f.Name), ...
                'FontSize', 7, 'Color', [0.5 0.5 0.5], 'Interpreter', 'none');
        end
    end
    title(sprintf('%d-Year BPT Time-Dependent Probability (3-Level Cascade Model)', forecastWindow));
    
    clim([0 maxPbpt]);
    cb = colorbar; ylabel(cb, 'BPT Conditional Probability');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
    % =====================================================================
    % MAP 5: PROBABILITY GAIN (Subsection Level 3D Patch Map)
    % =====================================================================
    figure('Name', 'NextQuake - Subsection Probability Gain Map', 'Color', 'w');
    hold on; grid on; axis equal; colormap(jet);
    
    for i = 1:numFaults
        f = faultResults(i);
        subIdx = f.GlobalSubIndices;
        if isempty(subIdx) || ~isfield(f, 'lon_mesh') || ~isfield(f, 'lat_mesh'), continue; end
        
        for idx = 1:length(subIdx)
            gID = subIdx(idx); 
            s_idx = subsectionData.S_idx(gID); 
            d_idx = subsectionData.D_idx(gID); 
            
            try
                lon_coords = [f.lon_mesh(s_idx, d_idx);     f.lon_mesh(s_idx+1, d_idx); ...
                              f.lon_mesh(s_idx+1, d_idx+1); f.lon_mesh(s_idx, d_idx+1)];
                          
                lat_coords = [f.lat_mesh(s_idx, d_idx);     f.lat_mesh(s_idx+1, d_idx); ...
                              f.lat_mesh(s_idx+1, d_idx+1); f.lat_mesh(s_idx, d_idx+1)];
                
                c = get_color_local(subGain(gID), subGain, jet);
                patch(lon_coords, lat_coords, c, 'EdgeColor', 'k', 'EdgeAlpha', 0.3, 'LineWidth', 0.4);
            catch
                continue;
            end
        end
        text(f.coordsWGS(1,1), f.coordsWGS(1,2), f.Name, 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');
    end
    title('Subsection Probability Gain Map (G = P_{BPT} / P_{Poisson})');
    cb = colorbar; ylabel(cb, 'Probability Gain Factor (G)');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    if max(subGain) > 0, clim([0 max(subGain)]); end

    % =====================================================================
    % MAP 6: PARTICIPATION RATE (Subsection Level 3D Patch Map - LOG10 SCALE)
    % =====================================================================
    figure('Name', 'NextQuake - Subsection Participation Rate Map (Log10)', 'Color', 'w');
    hold on; grid on; axis equal; colormap(jet);
    
    % Consideriamo solo i tassi strettamente positivi per evitare il log10(0)
    validIdx = subParticipationRate > 0;
    logSubParticipationRate = NaN(size(subParticipationRate));
    logSubParticipationRate(validIdx) = log10(subParticipationRate(validIdx));
    
    % Estraiamo i limiti reali dei logaritmi per la colorazione
    activeLogRates = logSubParticipationRate(validIdx);
    minLog = min(activeLogRates);
    maxLog = max(activeLogRates);
    
    for i = 1:numFaults
        f = faultResults(i);
        subIdx = f.GlobalSubIndices;
        if isempty(subIdx) || ~isfield(f, 'lon_mesh') || ~isfield(f, 'lat_mesh'), continue; end
        
        for idx = 1:length(subIdx)
            gID = subIdx(idx); 
            s_idx = subsectionData.S_idx(gID); 
            d_idx = subsectionData.D_idx(gID); 
            
            try
                lon_coords = [f.lon_mesh(s_idx, d_idx);     f.lon_mesh(s_idx+1, d_idx); ...
                              f.lon_mesh(s_idx+1, d_idx+1); f.lon_mesh(s_idx, d_idx+1)];
                          
                lat_coords = [f.lat_mesh(s_idx, d_idx);     f.lat_mesh(s_idx+1, d_idx); ...
                              f.lat_mesh(s_idx+1, d_idx+1); f.lat_mesh(s_idx, d_idx+1)];
                
                % Colorazione basata sul valore logaritmico effettivo
                if validIdx(gID)
                    % Usiamo la nuova funzione locale specifica per il logaritmo
                    c = get_color_local_log(logSubParticipationRate(gID), minLog, maxLog, jet);
                    patch(lon_coords, lat_coords, c, 'EdgeColor', 'k', 'EdgeAlpha', 0.3, 'LineWidth', 0.4);
                else
                    % Sotto-sezioni silenti (tasso = 0) colorate in grigio chiaro
                    patch(lon_coords, lat_coords, [0.9 0.9 0.9], 'EdgeColor', 'k', 'EdgeAlpha', 0.2, 'LineWidth', 0.4);
                end
            catch
                continue;
            end
        end
        text(f.coordsWGS(1,1), f.coordsWGS(1,2), f.Name, 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');
    end
    title('Total Subsection Event Participation Rates (Log10 Grid Projection)');
    cb = colorbar; ylabel(cb, 'Participation Rate [Log_{10}(events/yr)]');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
    % Force colorbar scale to match the log10 floors safely
    if ~isempty(activeLogRates) && maxLog > minLog
        clim([minLog maxLog]);
    end
end

% =====================================================================
% --- HELPER FUNCTIONS ---
% =====================================================================

function c = get_color_local(val, all_vals, cmap)
    % Helper function for linear maps (Maps 1 to 5)
    if isempty(all_vals) || max(all_vals) == min(all_vals) || val < 0 || isnan(val)
        c = cmap(1,:); return;
    end
    idx = round(1 + (size(cmap,1)-1) * (val - min(all_vals)) / (max(all_vals) - min(all_vals)));
    idx = max(1, min(size(cmap,1), idx));
    c = cmap(idx, :);
end

function c = get_color_local_log(val, min_val, max_val, cmap)
    % Helper function specific for negative/logarithmic scaling (Map 6)
    if isnan(val) || min_val == max_val
        c = cmap(1,:); return;
    end
    idx = round(1 + (size(cmap,1)-1) * (val - min_val) / (max_val - min_val));
    idx = max(1, min(size(cmap,1), idx));
    c = cmap(idx, :);
end
