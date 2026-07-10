function save_all_ruptures_SD_Batched(All_Final_Ruptures, subsectionData)
    % 1. Setup Output Directory
    output_dir = './Output/Ruptures';
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    
    total = height(All_Final_Ruptures);
    if total == 0, warning('No ruptures to save!'); return; end
    
    plots_per_fig = 12; % 3x4 grid
    num_figs = ceil(total / plots_per_fig);
    
    fprintf('💾 Saving %d batches (12 ruptures each) in %s...\n', num_figs, output_dir);
    
    for f = 1:num_figs
        % Create a NEW figure for each batch
        hFig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 1600 1200]);
        
        start_idx = (f-1) * plots_per_fig + 1;
        end_idx = min(f * plots_per_fig, total);
        
        % Create the 3x4 grid
        for i = start_idx:end_idx
            ax_idx = mod(i-1, plots_per_fig) + 1;
            subplot(3, 4, ax_idx);
            
            rupture_row = All_Final_Ruptures(i, :);
            
            plot_rupture_SD_view(rupture_row, subsectionData);
        end
        
        % Force rendering
        drawnow;
        
        % 2. Define Filename
        id_start = All_Final_Ruptures.Rupture_ID(start_idx);
        id_end   = All_Final_Ruptures.Rupture_ID(end_idx);
        filename = sprintf('Batch_Rup%d_to_Rup%d.png', id_start, id_end);
        full_path = fullfile(output_dir, filename);
        
        % 3. Save to PNG
        exportgraphics(hFig, full_path, 'Resolution', 150);
        
        close(hFig);
        fprintf('   ...saved batch %d/%d (%s)\n', f, num_figs, filename);
    end
    fprintf('✅ Done.\n');
end

%% ================== HELPER FUNCTION ==================
function plot_rupture_SD_view(rupture_row, Global_Subsections_Data)
    hold on; box on;
    
    % Identify involved elements
    rupture_IDs = rupture_row.Global_IDs{1};
    Rup_Sub = Global_Subsections_Data(ismember(Global_Subsections_Data.Global_ID, rupture_IDs), :);
    faults_involved = unique(Rup_Sub.Fault_ID);
    
    color_active = [0.85 0.33 0.1]; % orange
    color_bg = [0.92 0.92 0.92];    % light grey
    
    current_x_offset = 0;
    
    for f_id = faults_involved'
        f_data = Global_Subsections_Data(Global_Subsections_Data.Fault_ID == f_id, :);
        
        % Plot background
        scatter(f_data.S_idx + current_x_offset, -f_data.D_idx, 25, ...
            'MarkerFaceColor', color_bg, 'MarkerEdgeColor', [0.8 0.8 0.8], 'LineWidth', 0.3);
        
        % Plot active
        active_in_f = f_data(ismember(f_data.Global_ID, rupture_IDs), :);
        if ~isempty(active_in_f)
            scatter(active_in_f.S_idx + current_x_offset, -active_in_f.D_idx, 30, ...
                'MarkerFaceColor', color_active, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
        
        % Label Fault ID
        text(mean(f_data.S_idx) + current_x_offset, 0.5, sprintf('F%d', f_id), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');
            
        current_x_offset = current_x_offset + max(f_data.S_idx) + 2;
    end
    
    title(sprintf('Rup %d | Mw %.2f', rupture_row.Rupture_ID, rupture_row.Mw), 'FontSize', 9);
    set(gca, 'FontSize', 7, 'XTick', [], 'YTick', []);
    axis tight;
    yl = ylim; ylim([yl(1)-0.5, 1.5]);
end