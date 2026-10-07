function plotInversionAnalysisMaps(outputs, faultResults, subsectionData, invData, bestRates, params)
% PLOTINVERSIONANALYSISMAPS - Advanced spatial hazard visualizations optimized for 3D meshes:
% 1. Geological Target Mean Slip Rates
% 2. Mean Return Periods (Parent Fault Level)
% 3. Poisson Probabilities (Parent Fault Level)
% 4. BPT Time-Dependent Probabilities (Parent Fault Level)
% 5. Probability Gain Map (3D Subsection Patch Level)
% 6. Participation Rate Map (3D Subsection Patch Level in Log10 Scale)

    % --- Extract Parameters from Struct or Set Defaults ---
    if nargin < 6, params = struct(); end
    if isfield(params, 'forTwind'), forecastWindow = params.forTwind; else, forecastWindow = 30; end
    if isfield(params, 'forecastMw'), Mw_limit = params.forecastMw; else, Mw_limit = 6; end 
   
    numFaults = length(faultResults);
    numSubsections = height(subsectionData);
    
    % --- Extract Pre-Calculated Values from the Unified Output Engine ---
    faultParticipationRate = outputs.Faults.ParticipationRate;
    faultPoissonProb       = outputs.Faults.PoissonProb;
    faultBPTProb           = outputs.Faults.BPT_TimeDependentProb;
    
    subPoissonProb         = outputs.Subsections.PoissonProb;
    subBPTProb             = outputs.Subsections.BPT_TimeDependentProb;
    subParticipationRate   = outputs.Subsections.ParticipationRate;
    subGain                = outputs.Subsections.Gain;
    
    % --- Calculate Parent Fault Mean Return Periods ---
    returnPeriods = 1 ./ faultParticipationRate;
    returnPeriods(faultParticipationRate == 0) = NaN;
    
    fprintf('📊 Generating Spatial Analysis Maps (Mw >= %.1f, Forecast Window = %d years)...\n', Mw_limit, forecastWindow);    
    
    % =====================================================================
    % MAP 1: GEOLOGICAL TARGET SLIP RATES (Mean)
    % =====================================================================
    figure('Name', 'Target Slip Rates', 'Color', 'w');
    hold on; grid on; axis equal;
    colormap(flipud(hot));
    
    faultTargetMeanSR = zeros(numFaults, 1);
    if isfield(params, 'fileSR') && exist(params.fileSR, 'file')
        srTable = readtable(params.fileSR);
        for i = 1:numFaults
            fName = faultResults(i).Name;
            idx_table = find(strcmpi(srTable.FaultName, fName), 1);
            if ~isempty(idx_table)
                sMin = srTable.SRmin(idx_table);
                sMax = srTable.SRmax(idx_table); 
                faultTargetMeanSR(i) = (sMin + sMax) / 2;
            else
                warning('Fault %s not found in %s.', fName, params.fileSR);
                faultTargetMeanSR(i) = NaN;
            end
        end
    else
        warning('Input SR file not found or params.fileSR missing.');
    end
    
    validSR = faultTargetMeanSR(~isnan(faultTargetMeanSR));
    if ~isempty(validSR)
        minSR = min(validSR); maxSR = max(validSR);
    else
        minSR = 0; maxSR = 1;
    end
    
    for i = 1:numFaults
        f = faultResults(i);
        val = faultTargetMeanSR(i);
        if isnan(val), continue; end
        
        c = get_color_from_cmap(val, minSR, maxSR, flipud(hot));
        plot(f.coordsWGS(:,1), f.coordsWGS(:,2), '-', 'LineWidth', 3, 'Color', c);
        text(f.coordsWGS(1,1) + 0.01, f.coordsWGS(1,2), f.Name, 'FontSize', 7, 'Interpreter', 'none');
    end
    cb = colorbar; ylabel(cb, 'Target Mean Slip Rate [mm/yr]');
    clim([minSR, maxSR]); 
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
% =====================================================================
% MAP 2: RECURRENCE INTERVAL / RETURN PERIOD (Parent Fault Level)
% =====================================================================
figure('Name', 'Return Period', 'Color', 'w');
hold on; grid on; axis equal; colormap(jet);

validTR = returnPeriods(~isnan(returnPeriods) & returnPeriods < inf);
minTR = min(validTR); 
maxTR = max(validTR);

for i = 1:numFaults
    f = faultResults(i);
    tr_val = returnPeriods(i);
    
    if ~isnan(tr_val) && tr_val < inf
        % Plot line segment setting Z as the color value 
        % so MATLAB natively synchronizes the colorbar
        z_val = repmat(tr_val, size(f.coordsWGS(:,1)));
        patch([f.coordsWGS(:,1); NaN], [f.coordsWGS(:,2); NaN], [z_val; NaN], ...
              'EdgeColor', 'flat', 'LineWidth', 3, 'CData', [z_val; NaN]);
          
      
        text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf(' T_r: %.0f y', tr_val), ...
            'FontSize', 7, 'Interpreter', 'tex');
    else
        plot(f.coordsWGS(:,1), f.coordsWGS(:,2), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    end
end

title(sprintf('Mean Return Period (Years) for Mw \\geq %.1f', Mw_limit));
xlabel('Longitude (°)'); 
ylabel('Latitude (°)');

% Proper management of color limits and logarithmic scale
clim([minTR maxTR]);
set(gca, 'ColorScale', 'log');

cb = colorbar;
ylabel(cb, 'Return Period T_r (Years)', 'Interpreter', 'tex');

% =====================================================================
% MAP 3: POISSON PROBABILITY (Parent Fault Level)
% =====================================================================
maxP = max(faultPoissonProb(faultParticipationRate > 0));
if isempty(maxP) || maxP == 0, maxP = 1; end

figure('Name', sprintf('Poisson Probability %dy', forecastWindow), 'Color', 'w');
hold on; grid on; axis equal; colormap(parula);

for i = 1:numFaults
    f = faultResults(i);
    p_val = faultPoissonProb(i);
    
    if faultParticipationRate(i) > 0
        % Plot line segment setting Z as color data for automatic colorbar scaling
        z_val = repmat(p_val, size(f.coordsWGS(:,1)));
        patch([f.coordsWGS(:,1); NaN], [f.coordsWGS(:,2); NaN], [z_val; NaN], ...
              'EdgeColor', 'flat', 'LineWidth', 3, 'CData', [z_val; NaN]);
        
        text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf(' P_p: %.1f%%', p_val*100), ...
            'FontSize', 7, 'Interpreter', 'tex');
    else
        plot(f.coordsWGS(:,1), f.coordsWGS(:,2), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    end
end

title(sprintf('%d-Year Poisson Probability (Mw \\geq %.1f)', forecastWindow, Mw_limit));
xlabel('Longitude (°)'); 
ylabel('Latitude (°)');

% Proper color scale limits
clim([0 maxP]);

cb = colorbar; 
ylabel(cb, 'Probability');
    
    % =====================================================================
    % MAP 4: BPT TIME-DEPENDENT PROBABILITY (Parent Fault Level)
    % =====================================================================
    maxPbpt = max(faultBPTProb(faultParticipationRate > 0));
    if isempty(maxPbpt) || maxPbpt == 0, maxPbpt = 1; end
    
    figure('Name', sprintf('BPT Time-Dependent Probability %dy', forecastWindow), 'Color', 'w');
    hold on; grid on; axis equal; colormap(spring);
    
    for i = 1:numFaults
        f = faultResults(i);
        if faultParticipationRate(i) > 0
           
            c = get_color_from_cmap(faultBPTProb(i), 0, maxPbpt, spring);
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), '-', 'LineWidth', 3, 'Color', c);
            text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf('P_{BPT}: %.1f%%', faultBPTProb(i)*100), ...
                'FontSize', 7, 'Interpreter', 'tex');
        else
            plot(f.coordsWGS(:,1), f.coordsWGS(:,2), ':', 'Color', [0.7 0.7 0.7]);
            text(f.coordsWGS(1,1), f.coordsWGS(1,2), sprintf('P_{BPT}: 0.0%%'), ...
                'FontSize', 7, 'Color', [0.5 0.5 0.5], 'Interpreter', 'tex');
        end
    end
    title(sprintf('BPT Time-Dependent Probability in the next %d-Year', forecastWindow));
    clim([0 maxPbpt]);
    cb = colorbar; ylabel(cb, 'Probability');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
    % =====================================================================
    % MAP 5: PROBABILITY GAIN (Subsection Level 3D Patch Map)
    % =====================================================================
    figure('Name', 'Subsection Probability Gain Map', 'Color', 'w');
    hold on; grid on; axis equal; colormap(jet);
    
    minGain = min(subGain);
    maxGain = max(subGain);
    if minGain == maxGain, maxGain = minGain + 1; end
    
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
                
                c = get_color_from_cmap(subGain(gID), minGain, maxGain, jet);
                patch(lon_coords, lat_coords, c, 'EdgeColor', 'k', 'EdgeAlpha', 0.3, 'LineWidth', 0.4);
            catch
                continue;
            end
        end
        % text(f.coordsWGS(1,1), f.coordsWGS(1,2), f.Name, 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');
    end
    title('Subsection Probability Gain Map (G = P_{BPT} / P_{Poisson})');
    cb = colorbar; ylabel(cb, 'Probability Gain (G)');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    clim([minGain maxGain]);
    
    % =====================================================================
    % MAP 6: PARTICIPATION RATE (Subsection Level 3D Patch Map - LOG10 SCALE)
    % =====================================================================
    figure('Name', 'Subsection Participation Rate Map', 'Color', 'w');
    hold on; grid on; axis equal; colormap(jet);
    
    validIdx = subParticipationRate > 0;
    logSubParticipationRate = NaN(size(subParticipationRate));
    logSubParticipationRate(validIdx) = log10(subParticipationRate(validIdx));
    
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
                
                if validIdx(gID)
                    c = get_color_local_log(logSubParticipationRate(gID), minLog, maxLog, jet);
                    patch(lon_coords, lat_coords, c, 'EdgeColor', 'k', 'EdgeAlpha', 0.3, 'LineWidth', 0.4);
                else
                    patch(lon_coords, lat_coords, [0.9 0.9 0.9], 'EdgeColor', 'k', 'EdgeAlpha', 0.2, 'LineWidth', 0.4);
                end
            catch
                continue;
            end
        end
        % text(f.coordsWGS(1,1), f.coordsWGS(1,2), f.Name, 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');
    end
    title('Total Subsection Event Participation Rates');
    cb = colorbar; ylabel(cb, 'Participation Rate [Log_{10}(events/yr)]');
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    
    if ~isempty(activeLogRates) && maxLog > minLog
        clim([minLog maxLog]);
    end
end

% =====================================================================
% --- HELPER FUNCTIONS ---
% =====================================================================
function c = get_color_local_log(val, min_val, max_val, cmap)
    if isnan(val) || min_val == max_val
        c = cmap(1,:); return;
    end
    idx = round(1 + (size(cmap,1)-1) * (val - min_val) / (max_val - min_val));
    idx = max(1, min(size(cmap,1), idx));
    c = cmap(idx, :);
end

function c = get_color_from_cmap(val, minVal, maxVal, cmap)
    if maxVal == minVal || isnan(val)
        idx = 1;
    else
        normVal = (val - minVal) / (maxVal - minVal);
        normVal = max(0, min(1, normVal));
        idx = round(normVal * (size(cmap, 1) - 1)) + 1;
    end
    c = cmap(idx, :);
end
