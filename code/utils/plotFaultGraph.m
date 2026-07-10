function plotFaultGraph(faultResults, subsectionData, A_S2S)
    % PLOTFAULTGRAPH_DARK - Generates a high-visibility connectivity graph
    % Purpose: Visualizes how 3D subsections connect individual faults into systems.
    
    %% 1. DATA EXTRACTION & MAPPING
    fIDs = [faultResults.Fault_ID];
    numFaults = numel(fIDs);
    Adj_Faults = zeros(numFaults); % Fault-to-fault adjacency matrix
    
    % Map Subsection-to-Subsection connections to Fault-to-Fault connections
    [subIdxA, subIdxB] = find(A_S2S);
    for k = 1:length(subIdxA)
        idA = subsectionData.Fault_ID(subIdxA(k));
        idB = subsectionData.Fault_ID(subIdxB(k));
        
        idxA = find(fIDs == idA);
        idxB = find(fIDs == idB);
        
        if ~isempty(idxA) && ~isempty(idxB) && idxA ~= idxB
            Adj_Faults(idxA, idxB) = 1;
            Adj_Faults(idxB, idxA) = 1;
        end
    end
    
    % Create the Graph object
    G = graph(Adj_Faults);
    
    %% 2. FIGURE SETUP (DARK MODE)
    bg_color = [0.05 0.05 0.05]; % Very dark grey for high contrast
    fig = figure('Name', 'NEXTQUAKE: Fault Network Topology', 'Color', bg_color, ...
                 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    
    % Initialize axes with matching background
    ax = axes('Parent', fig, 'Color', bg_color, 'XColor', 'none', 'YColor', 'none');
    hold(ax, 'on');

    %% 3. GRAPH PLOTTING
    % Note: We plot WITHOUT labels first to handle centering manually
    p = plot(G, 'Layout', 'force', ...
               'WeightEffect', 'direct', ...
               'Iterations', 1000, ...
               'NodeLabel', {}); % Suppress default labels

    % --- Node Styling (Electric Blue) ---
    p.Marker = 'o';
    p.MarkerSize = 28;           % Large nodes to accommodate numbers
    p.NodeColor = [0.0 0.5 1.0]; % Neon blue
    
    % --- Edge Styling (Subtle Grey) ---
    p.EdgeColor = 'w'; % Light grey for visibility on black
    p.LineWidth = 1.5;
    p.EdgeAlpha = 0.4;           % Transparency to reduce visual clutter
    
    %% 4. SPATIAL EXPANSION & LABEL CENTERING
    % Increase node spacing to avoid overlaps
    p.XData = p.XData * 3.0; 
    p.YData = p.YData * 3.0;
    
    % Manually place ID numbers at the center of each node
    for i = 1:numFaults
        text(p.XData(i), p.YData(i), num2str(fIDs(i)), ...
            'Color', 'w', ...            % White text
            'FontSize', 10, ...          % Adjust size for 2-digit IDs
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Consolas');     % Monospaced font for better centering
    end

    %% 5. ANNOTATIONS & FORMATTING
    total_conns = sum(Adj_Faults(:))/2;
    title_str = sprintf('FAULT CONNECTIVITY NETWORK\nNodes: %d | Inter-Fault Bridges: %d', ...
                        numFaults, total_conns);
    
    title(ax, title_str, 'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Optimization for display
    axis(ax, 'equal'); 
    axis(ax, 'off');
    set(ax, 'LooseInset', get(ax, 'TightInset'));
    
    hold(ax, 'off');
    fprintf('✅ Dark Mode Graph generated with %d connections.\n', total_conns);
end