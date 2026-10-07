function plotS2SConnections(faultResults, Global_Subsections_Data, A_S2S)
% PLOTS2SCONNECTIONS - Visualizes Subsection-to-Subsection (S2S) propagation bridges.
%
% This plot helps validate the connectivity matrix (A_S2S). It shows
% nodes (centroids) connected by lines wherever a rupture jump is possible.

    fprintf('📊 Plotting S2S connections...\n');

    if isempty(Global_Subsections_Data) || isempty(A_S2S)
        warning('Missing data or S2S matrix. Cannot plot.');
        return;
    end

    figure('Name', 'S2S Propagation Bridges', 'Color', 'w');
%     set(gca, 'Color', [0.1 0.1 0.1]); % Sfondo grigio scuro/nero
% set(gcf, 'Color', [0.1 0.1 0.1]);
% grid off;
    hold on; grid on; axis equal;

    % --- Step 1: Prepare Colors ---
    % Group subsections by Fault_ID to assign a distinct color per fault
    faultIDs = Global_Subsections_Data.Fault_ID;
    uniqueFIDs = unique(faultIDs);
    nFaults = length(uniqueFIDs);
    cmap = colorcube(nFaults); % Uses colorcube for very distinct colors

    % --- Step 2: Plot S2S Connections (Black lines) ---
    % Find all pairs (i, j) where A_S2S is 1
    [rows, cols] = find(A_S2S);
    
    % Get Centroids (X, Y)
    subX = Global_Subsections_Data.X_km;
    subY = Global_Subsections_Data.Y_km;

    % Iterate over connections (only need the upper triangle to avoid duplicates)
    for k = 1:length(rows)
        i = rows(k);
        j = cols(k);
        
        if i < j % Plot each connection once
            % Draw a black line connecting centroid i to centroid j
            plot([subX(i), subX(j)], [subY(i), subY(j)], 'k-', ...
                'LineWidth', 1.2, 'HandleVisibility', 'off');
        end
    end

    % --- Step 3: Plot Nodes (Centroids) on top ---
    % We plot nodes *after* lines to make them appear on top
    node_handles = [];
    fault_names = {};

    for f = 1:nFaults
        fid = uniqueFIDs(f);
        idx = find(faultIDs == fid);
        
        if isempty(idx), continue; end
        
        % Plot nodes (centroids)
        h = plot(subX(idx), subY(idx), 'o', 'MarkerFaceColor', cmap(f,:), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 5);
        
        % Collect handles for legend (limit to a reasonable number)
        if f <= 100 
            node_handles(end+1) = h;
            fault_names{end+1} = faultResults(fid).Name;
        end
    end

    % --- Step 4: Formatting ---
    xlabel('Easting (UTM km)');
    ylabel('Northing (UTM km)');
    title('S2S nearest-only propagation bridges (Only distinct faults)');
    
    if ~isempty(node_handles)
        legend(node_handles, fault_names, 'Location', 'northeastoutside', 'Interpreter', 'none');
    end
    
    hold off;
end
