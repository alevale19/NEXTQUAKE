function plotCentroids3D(faultResults, subsectionData)
% PLOTCENTROIDS3D - Plots 3D centroids and surface fault traces
%
% Inputs:
%    faultResults   - Struct with original traces (coordsWGS)
%    subsectionData - Table with X_km, Y_km, Z_km

    figure('Name', '3D Fault System', 'Color', 'w');
    hold on;

    % 1. Plot Surface Traces (Blue lines at Z = 0)
    % Since the table uses UTM km, we must convert WGS84 traces to UTM km
    for nf = 1:length(faultResults)
        f = faultResults(nf);
        
        % Convert WGS84 to UTM (meters)
        [E, N, ~] = deg2utm(f.coordsWGS(:,2), f.coordsWGS(:,1));
        
        % Plot in Km at Z=0
        plot3(E/1000, N/1000, zeros(size(E)), 'b-', 'LineWidth', 1.5);
    end

    % 2. Plot Subsections Centroids
    % Using scatter3: X, Y, -Z (depth)
    % Color is mapped to depth (Z_km)
    s1 = scatter3(subsectionData.X_km, subsectionData.Y_km, -subsectionData.Z_km, ...
                  25, subsectionData.Z_km, 'filled', 'MarkerEdgeColor', 'k');

    % 3. Formatting
    colorbar;
    ylabel(colorbar, 'Depth (km)');
    colormap(flipud(hot)); % Hot colors for deeper sections
    
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    zlabel('Depth (km)');
    title('3D View: Fault Surface Traces & Subsection Centroids');
    
    axis equal;
    grid on;
    view(3); % Standard 3D perspective view
    
    hold off;
end
