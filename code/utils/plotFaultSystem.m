function plotFaultSystem(faultResults, params)
% PLOTFAULTSYSTEM - Generates two figures to visualize fault geometry in a 2D map.
% This version keeps the 2D view and correctly shows the mesh box starting 
% from the upper seismogenic edge.

    if isempty(faultResults)
        warning('No fault data available to plot.');
        return;
    end
    
    % =========================================================
    % FIGURE 1: 2D MAP GEOMETRY & MESH BOX (FROM UPPER EDGE)
    % =========================================================
    figure('Name', 'Fault Discretization Map', 'Color', 'w');
    hold on; grid on; axis equal;
    
    h1 = []; h2 = []; h3 = []; h4 = []; % Handles for structural legend
    for nf = 1:length(faultResults)
        f = faultResults(nf);
        
        % 1. Plot Main Surface Trace (Blue - real surface alignment Z=0)
        h1 = plot(f.coordsWGS(:,1), f.coordsWGS(:,2), 'b-', 'LineWidth', 2);
        
        % 2. Plot 2D Mesh Box (Discretization grid starting from Upper Edge)
        if isfield(f, 'lat_mesh') && isfield(f, 'lon_mesh')
          
            h4 = plot(f.lon_mesh, f.lat_mesh, 'k-', 'LineWidth', 0.5);
            plot(f.lon_mesh', f.lat_mesh', 'k-', 'LineWidth', 0.5);    
         
            % --- Bottom Trace ---
            h2 = plot(f.lon_mesh(:, end), f.lat_mesh(:, end), 'r--', 'LineWidth', 1);
        end
        
        % 3. Plot Upper Edge Nodes (Yellow markers on the start of the mesh box)
        
        if isfield(f, 'lon_mesh') && isfield(f, 'lat_mesh')
            h3 = plot(f.lon_mesh(:, 1), f.lat_mesh(:, 1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 4);
        end
        
        % Fault Name Label
        text(f.coordsWGS(1,1), f.coordsWGS(1,2), f.Name, 'FontSize', 8, ...
            'FontWeight', 'bold', 'Interpreter', 'none');
    end
    
    % --- Overlay Constraint Data (Field Sites) ---
    h_paleo = []; h_point = [];
    
    if isfield(params, 'filePaleo') && exist(params.filePaleo, 'file')
        pal = readtable(params.filePaleo);
        if any(ismember(pal.Properties.VariableNames, {'Lon', 'Lat'}))
            h_paleo = plot(pal.Lon, pal.Lat, 'gd', 'MarkerSize', 8, ...
                'MarkerFaceColor', 'g', 'LineWidth', 1.2);
        end
    end
    
    if isfield(params, 'filePointSR') && exist(params.filePointSR, 'file')
        pSR = readtable(params.filePointSR);
        if any(ismember(pSR.Properties.VariableNames, {'Lon', 'Lat'}))
            h_point = plot(pSR.Lon, pSR.Lat, 'ms', 'MarkerSize', 7, ...
                'MarkerFaceColor', 'm', 'LineWidth', 1.2);
        end
    end
    
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    title('Fault System Geometry');
    
    % Safe Legend Management for Figure 1
    f1_handles = []; f1_labels = {};
    if ~isempty(h1), f1_handles(end+1) = h1; f1_labels{end+1} = 'Surface Trace (Z=0)'; end
    if ~isempty(h2), f1_handles(end+1) = h2; f1_labels{end+1} = 'Bottom Trace (lowerST)'; end
    if ~isempty(h3), f1_handles(end+1) = h3; f1_labels{end+1} = 'Upper Edge Nodes (upperST)'; end
    if ~isempty(h4), f1_handles(end+1) = h4(1); f1_labels{end+1} = 'Mesh Box (Subsections)'; end
    if ~isempty(h_paleo), f1_handles(end+1) = h_paleo; f1_labels{end+1} = 'Paleo Sites'; end
    if ~isempty(h_point), f1_handles(end+1) = h_point; f1_labels{end+1} = 'Point Slip Rates'; end
    
    if ~isempty(f1_handles)
        legend(f1_handles, f1_labels, 'Location', 'northeastoutside'); 
    end
    hold off;
    
    % =========================================================
    % FIGURE 2: SYSTEM CONNECTIVITY (FAUL CLUSTERS)
    % =========================================================
    figure('Name', 'Fault Clusters', 'Color', 'w');
    hold on; grid on; axis equal; 
    
    allIDs = [faultResults.SystemID];
    uIDs = unique(allIDs);
    nSystems = length(uIDs);
    cmap = lines(nSystems); 
    
    h_sys = [];
    for nf = 1:length(faultResults)
        f = faultResults(nf);
        colorIdx = find(uIDs == f.SystemID);
        sysColor = cmap(colorIdx, :);

        h_sys = plot(f.coordsWGS(:,1), f.coordsWGS(:,2), '-', 'Color', sysColor, ...
            'LineWidth', 3);
    end
    
    xlabel('Longitude (°)'); ylabel('Latitude (°)');
    title(sprintf('Fault Clusters: %d Independent Systems Identified', nSystems));
    
    % if ~isempty(h_sys)
    %     legend(h_sys(1), 'System Cluster', 'Location', 'northeastoutside');
    % end
    hold off;
end
