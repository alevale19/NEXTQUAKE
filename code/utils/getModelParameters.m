function params = getModelParameters()
% GETMODELPARAMETERS - Aligned 3-column layout.
    import javax.swing.*
    import java.awt.*
    
    % --- Main Dialog Panel ---
    mainPanel = JPanel(GridBagLayout());
    gbc = GridBagConstraints();
    gbc.insets = Insets(10, 15, 10, 15); 
    gbc.anchor = GridBagConstraints.NORTHWEST;
    gbc.fill = GridBagConstraints.HORIZONTAL;
    
    % Fonts
    titleFont = Font('Arial', Font.BOLD, 18);
    labelFont = Font('Arial', Font.PLAIN, 18);
    
    % =====================================================================
    % COLUMN 1: RUPTURE CATALOG SETUP
    % =====================================================================
    col1 = JPanel();
    col1.setLayout(BoxLayout(col1, BoxLayout.Y_AXIS));
    
    t1 = JLabel('<html>RUPTURE CATALOG SETUP</html>'); 
    t1.setFont(titleFont); t1.setForeground(Color(0, 0.4, 0.7));
    col1.add(t1); col1.add(Box.createRigidArea(Dimension(0, 10)));
    
    labels1 = {'Input Faults file:', 'Forecast Mw:', 'Fault Mechanism:', 'Regional b-value:', 'Max Gap Shallow (km):', ...
               'Max Gap Deep (km):', 'Z Lim (km):', 'Max Azimuth Diff (°):', ...
               'Max Rake Diff (°):', 'Max Dip Diff (°):', 'Max Aspect Ratio (AR):', ...
               'Max Mw:', 'Output Catalog Name:'};
    defs1   = {'InputFaults', '5.9', 'N', '1.0', '6.0', '6.0', '5.0', '40', '30', '20', '8.0', '7.7', 'RupCatalog'};
    
    fields1 = cell(1, length(labels1));
    for i = 1:length(labels1)
        lbl = JLabel(labels1{i}); lbl.setFont(labelFont);
        col1.add(lbl);
        fields1{i} = JTextField(defs1{i}, 18);
        fields1{i}.setMaximumSize(Dimension(300, 35));
        fields1{i}.setAlignmentX(Component.LEFT_ALIGNMENT);
        col1.add(fields1{i});
        col1.add(Box.createRigidArea(Dimension(0, 6)));
    end
    
    % =====================================================================
    % COLUMN 2: INVERSION DATA & CONSTRAINTS (Updated for Time Settings)
    % =====================================================================
    col2 = JPanel();
    col2.setLayout(BoxLayout(col2, BoxLayout.Y_AXIS));
    
    t2a = JLabel('INVERSION DATA & CONSTRAINTS');
    t2a.setFont(titleFont); t2a.setForeground(Color(0.8, 0.3, 0));
    col2.add(t2a); col2.add(Box.createRigidArea(Dimension(0, 10)));
    
    % Added 'Fault Time Settings (Opt):' at index 8
    labels2a = {'Slip Rate:', 'Paleoseismic Rate (Opt):', 'Regional MFD (Opt):', ...
                'Fault MFD (Opt):', 'Segmentation Constraint (Opt):', 'Point Slip Rate (Opt):', ...
                'Pre-processed Rupture Catalog (.mat):', 'Fault Time Settings (Opt):', 'Forecast Time Window:','Year for computing elapsed time'};
    defs2a   = {'SR', '', '', '', '', '', '', 'Fault_TimeSettings', '50','2026'};
    
    fields2a = cell(1, length(labels2a));
    for i = 1:length(labels2a)
        lbl = JLabel(labels2a{i}); lbl.setFont(labelFont);
        col2.add(lbl);
        fields2a{i} = JTextField(defs2a{i}, 18);
        fields2a{i}.setMaximumSize(Dimension(300, 35));
        fields2a{i}.setAlignmentX(Component.LEFT_ALIGNMENT);
        col2.add(fields2a{i});
        col2.add(Box.createRigidArea(Dimension(0, 6)));
    end
    
    % =====================================================================
    % COLUMN 3: SIMULATED ANNEALING SETTINGS
    % =====================================================================
    col3 = JPanel();
    col3.setLayout(BoxLayout(col3, BoxLayout.Y_AXIS));
    
    t2b = JLabel('SIMULATED ANNEALING SETTINGS');
    t2b.setFont(titleFont); t2b.setForeground(Color(0.2, 0.5, 0.2));
    col3.add(t2b); col3.add(Box.createRigidArea(Dimension(0, 10)));
    
    labels2b = {'Weight Paleoseismic Rate:', 'Weight Regional MFD:', ...
                'Weight Fault MFD:', 'Weight Segmentation:', 'SA Runs:', ...
                'Number of Iterations:'};
    defs2b   = {'', '0', '', '', '2', '1000000'};
    
    fields2b = cell(1, length(labels2b));
    for i = 1:length(labels2b)
        lbl = JLabel(labels2b{i}); lbl.setFont(labelFont);
        col3.add(lbl);
        fields2b{i} = JTextField(defs2b{i}, 18);
        fields2b{i}.setMaximumSize(Dimension(300, 35));
        fields2b{i}.setAlignmentX(Component.LEFT_ALIGNMENT);
        col3.add(fields2b{i});
        col3.add(Box.createRigidArea(Dimension(0, 6)));
    end
    
    % --- Add Columns to GridBagLayout ---
    gbc.gridy = 0;
    gbc.gridx = 0; mainPanel.add(col1, gbc);
    gbc.gridx = 1; mainPanel.add(col2, gbc);
    gbc.gridx = 2; mainPanel.add(col3, gbc);
    
    % --- Display Buttons Dialog ---
    options = {'RUN NEXTQUAKE', 'CANCEL'};
    choice = JOptionPane.showOptionDialog([], mainPanel, ...
        'NEXTQUAKE Setup', ...
        JOptionPane.YES_NO_OPTION, ...
        JOptionPane.PLAIN_MESSAGE, ...
        [], options, options{1});
        
    if choice ~= 0
        fprintf('⚠️ Execution canceled by user.\n');
        params = []; return;
    end
    
    % --- DATA PARSING ---
    try
        % Extract from Column 1
        params.forecastMw        = str2double(char(fields1{2}.getText()));
        params.faultMechanism    = strtrim(char(fields1{3}.getText())); % N, SS, R, All, SUB
        params.b_val             = str2double(char(fields1{4}.getText()));
        params.maxGapShallow     = str2double(char(fields1{5}.getText()));
        params.maxGapDeep        = str2double(char(fields1{6}.getText()));
        params.zLimit            = str2double(char(fields1{7}.getText()));
        params.maxAzDiff         = str2double(char(fields1{8}.getText()));
        params.maxRakeDiff       = str2double(char(fields1{9}.getText()));
        params.maxDipDiff        = str2double(char(fields1{10}.getText()));
        params.maxAspectRatio    = str2double(char(fields1{11}.getText()));
        params.maxMw             = str2double(char(fields1{12}.getText()));
        
        % other parameters not yet in the GUI
        params.coordSystem = 'WGS84';
        params.maxGap = max(params.maxGapShallow, params.maxGapDeep);
        params.minSubs = 1;
        params.applyProbVisible = 1; % 1 (yes) or 0 (no). It applies Pizza et al., 2023 eqs to Gsr matrix
        params.mu = 3e10; % Shear Modulus in Pascal (Nm^-2)
        params.MoRateReduction = 0.0;  % Optional additional moment rate reduction (0 = none, 1 = full reduction). It applies only to regional MFD.
        params.mag_delta = 0.1;        % Magnitude bin width for building the truncated GR used as target.
        params.weightL2          = 0.0001; % Weight Smoothing (L2)
        
        %params.seed = randi(100000); %uncomment to get a seed number of each run
        
        params.minRecommendedK = 500; % average number of times that the SA explores a single rupture

        % Extract from Column 2 (Note the index shift: 8 is now Time Settings, 9 is Window)
        params.forTwind          = str2double(char(fields2a{9}.getText()));
        params.currentYear = str2double(char(fields2a{10}.getText()));

        % Extract from Column 3
        params.weightPaleo       = str2double(char(fields2b{1}.getText()));
        params.weightMFD         = str2double(char(fields2b{2}.getText())); 
        params.weightLocalMFD    = str2double(char(fields2b{3}.getText())); 
        params.weightSeg         = str2double(char(fields2b{4}.getText()));
        params.runSA             = str2double(char(fields2b{5}.getText()));
        params.numIter           = str2double(char(fields2b{6}.getText()));
        
        % Paths and Directories Setup
        inputDir = './Input/';
        constrDir = fullfile(inputDir, 'Constraints');
        
        params.Faultdata = readtable(fullfile(inputDir, [char(fields1{1}.getText()), '.txt']));
        params.fileSR = fullfile(constrDir, [char(fields2a{1}.getText()), '.txt']);
        params.saveFile = [char(fields1{13}.getText()), '.txt'];
        
        % --- NEW: Parse and Load Time Settings Table ---
        timeFileName = char(fields2a{8}.getText());
        timeFilePath = fullfile(inputDir, [timeFileName, '.txt']);
        
        if ~isempty(timeFileName) && exist(timeFilePath, 'file')
            params.timeSettings = readtable(timeFilePath);
            fprintf('📅 Time-dependent settings file loaded successfully.\n');
        else
            params.timeSettings = table(); % Empty table fallback if file doesn't exist or is omitted
            fprintf('ℹ️ No time-dependent settings file found or specified. Fallback to default BPT/Poisson.\n');
        end
        
        % Optional Pre-processed .mat Catalog
        preloadedName = char(fields2a{7}.getText());
        if isempty(preloadedName), params.preloadedData = ''; else, params.preloadedData = fullfile(inputDir, [preloadedName, '.mat']); end
       
        % Optional Constraints Files mapping
        paleoName = char(fields2a{2}.getText());
        if isempty(paleoName), params.filePaleo = ''; else, params.filePaleo = fullfile(constrDir, [paleoName, '.txt']); end
        
        mfdName = char(fields2a{3}.getText()); 
        if isempty(mfdName), params.fileMFD = ''; else, params.fileMFD = fullfile(constrDir, [mfdName, '.txt']); end
        
        LocalmfdName = char(fields2a{4}.getText()); 
        if isempty(LocalmfdName), params.fileLocalMFD = ''; else, params.fileLocalMFD = fullfile(constrDir, [LocalmfdName, '.txt']); end
        
        segConName = char(fields2a{5}.getText());
        if isempty(segConName), params.fileSeg = ''; else, params.fileSeg = fullfile(constrDir, [segConName, '.txt']); end
        
        srPointName = char(fields2a{6}.getText());
        if isempty(srPointName), params.filePointSR = ''; else, params.filePointSR = fullfile(constrDir, [srPointName, '.txt']); end
        
        fprintf('🚀 Parameters set. Ready to compute.\n');
    catch ME
        fprintf('❌ Error parsing parameters: %s\n', ME.message);
        params = [];
    end
end
