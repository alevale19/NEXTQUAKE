%% NEXTQUAKE - Main Processing Script
clc; clear; close all;
set(groot, 'DefaultFigureWindowStyle', 'docked');
addpath('./utils');

logFileName = sprintf('NEXTQUAKE_Log_%s.txt', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
diary(logFileName); % Avvia la registrazione di tutto l'output della Command Window
diary on;

%% 1. INITIALIZATION & PARAMETERS
params = getModelParameters();
if isempty(params), return; end

% Define common paths for saving/loading
[filepath, name, ~] = fileparts(params.saveFile);
invDataPath = fullfile(filepath, [name, '_InvData.mat']);

%% 2. DATA SOURCE DECISION (Load vs Generate)
if ~isempty(params.preloadedData)
    fprintf('📂 Loading specified pre-processed data: %s\n', params.preloadedData);
    load(params.preloadedData);
    runFullPipeline = false;
else
    runFullPipeline = true;
end

%% 3. RUPTURE GENERATION PIPELINE
if runFullPipeline
    fprintf('🚀 Starting FULL NEXTQUAKE processing pipeline...\n');
    
    % --- A. Geometry & Discretization ---
    L_sub = seismicScalingEngine(params.forecastMw, 'Mw2Length', params.faultMechanism);
    faultResults = processFaultDatabase(params.Faultdata, params, L_sub);
    if isempty(faultResults), error('❌ Fault processing failed.'); end
    subsectionData = computeSubsectionCentroids(faultResults, params.timeSettings, params);
    
    % --- B. Connectivity & Graph Building ---
    faultResults = groupFaultSystems(faultResults, params);
    A_S2S = computeS2SMatrix(faultResults, subsectionData, params);
    
    % --- C. Rupture Generation ---
    fprintf('\n🧱 Generating Rupture Catalog...\n');
    genTimer = tic;
    All_Elementary_Ruptures = generateElementaryRuptures(faultResults, subsectionData, params);
    All_Final_Ruptures = generateMultiFaultRuptures(All_Elementary_Ruptures, subsectionData, A_S2S, params);
    fprintf('✅ Catalog Generation Complete in %.2f seconds.\n', toc(genTimer));
    
    % --- D. Prepare Data for Inversion ---
    cols = {'Mw', 'Area', 'Global_IDs'};
    if isempty(All_Final_Ruptures) || height(All_Final_Ruptures) == 0
        fprintf('ℹ️ No complex multi-fault ruptures found. Using elementary ruptures only.\n');
        All_Rupture_Candidates = All_Elementary_Ruptures(:, cols);
    else
        All_Rupture_Candidates = [All_Elementary_Ruptures(:, cols); All_Final_Ruptures(:, cols)];
    end
    
    All_Rupture_Candidates.Rup_Index = (1:height(All_Rupture_Candidates))';
    
    % Map subsections and calculate total fault area early
    for i = 1:length(faultResults)
        matchedRows = strcmp(subsectionData.FaultName, faultResults(i).Name);
        faultResults(i).GlobalSubIndices = subsectionData.Global_ID(matchedRows);
        faultResults(i).TotalArea_km2 = sum(subsectionData.Area_km2(matchedRows));
    end
    
    % Load Constraints (Slip Rates, Paleo, MFD from files...)
    constraints = loadInversionConstraints(subsectionData, params);
    invData = prepareInversionData(All_Rupture_Candidates, subsectionData, constraints, params);
    
    % --- E. Save for Future Use ---
    save(invDataPath, 'invData', 'subsectionData', 'faultResults', 'params', 'A_S2S', 'All_Rupture_Candidates', 'constraints', '-v7.3');
    fprintf('💾 Pre-processed data saved to: %s\n', invDataPath);
end

%% 4. INVERSION CONSTRAINTS & TARGETS
% Setup Inversion Targets
targets = setupInversionTargets(invData, subsectionData, All_Rupture_Candidates, constraints, params);

% Sanitize MFD Targets (Prevent NaN/out-of-bounds propagation)
if isfield(targets, 'mfdTarget') && any(isnan(targets.mfdTarget))
    warning('⚠️ NaN detected in targets.mfdTarget. Replacing NaN entries with 0.');
    targets.mfdTarget(isnan(targets.mfdTarget)) = 0;
end

%% 5. INITIALIZE RATES & TEST FAST ENERGY INTEGRITY
rates = zeros(invData.numRuptures, 1);

% Initialize state for Fast Delta-Energy Updating engine
initialState = initEnergyState(rates, invData, targets);

fprintf('\n⚡ Initial Energy Test:\n');
fprintf('   Total Energy: %.2e (SR: %.2e, MFD: %.2e, Paleo: %.2e, Seg: %.2e, L2: %.2e)\n\n', ...
        initialState.totalEnergy, initialState.energySR, initialState.energyMFD, ...
        initialState.energyPaleo, initialState.energySeg, initialState.energyL2);

if isnan(initialState.totalEnergy)
    error('❌ Initial Energy computed as NaN! Please verify targets.mfdTarget and rupture Mw boundaries.');
end

%% 6. LAUNCH SIMULATED ANNEALING
[bestRates, energyHistory] = runSimulatedAnnealing(invData, targets, params, faultResults);

%% 7. QUICK RESULTS CHECK
% Fig. 4 Convergence History
figure('Name', 'Convergence History', 'Color', 'w');
semilogy(energyHistory, 'LineWidth', 1.5);
grid on;
xlabel('Iteration'); ylabel('Total Energy');
title('SA Convergence');

% Fig. 5 Event Rate Match (Paleoseismic Comparison)
finalEventRates = sum(invData.Gsr' .* bestRates.mean', 2);
figure('Name', 'Event Rate Match', 'Color', 'w');
hold on;

if isfield(invData, 'paleoTargetIdx') && ~isempty(invData.paleoTargetIdx)
    pIdx = invData.paleoTargetIdx;
    pRate = targets.paleoRate;
    pStd = targets.paleoStd;
    
    errorbar(pIdx, pRate, pStd, 'go', 'MarkerFaceColor', 'g', ...
        'LineWidth', 1.2, 'DisplayName', 'Paleo Targets (1\sigma)');
    plot(pIdx, finalEventRates(pIdx), 'rs', 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Model at Paleo Sites');
end

plot(finalEventRates, 'r-', 'LineWidth', 1, 'DisplayName', 'Model (mean)');
legend('show', 'Location', 'northeast');
grid on; box on;
title('Event Rate Inversion Results & Paleoseismic Match');
ylabel('Event Rate (events/yr)');
xlabel('Subsection Index');
set(gca, 'YScale', 'log');

%% 8. VISUALIZATION & QUALITY CHECK
fprintf('\n📊 Launching Visualizations...\n');

% Geometry, Clusters & Connectivity
plotFaultSystem(faultResults, params); 
plotCentroids3D(faultResults, subsectionData); 
plotS2SConnections(faultResults, subsectionData, A_S2S); 

% Magnitude-Frequency Distribution
plotMFDComparison(bestRates, invData, targets, params); 

% Tectonic Moment Budget Check
checkMomentBalance(bestRates, invData, subsectionData, targets, params); 

% Post-processing BPT & Hazard Output Maps
outputs = processAnnealingResults(bestRates, invData, faultResults, subsectionData, params);
plotInversionAnalysisMaps(outputs, faultResults, subsectionData, invData, bestRates, params);

% Multi-Fault Cascades & Interaction Diagnostics
analyzeFaultCoRupturePairs(faultResults, subsectionData, invData, bestRates);
analyzeComplexMultiFaultClusters(faultResults, subsectionData, invData, bestRates);
reportMultiFaultContribution(faultResults, subsectionData, invData, bestRates);

% Global Participation Map
plotGlobalParticipationM_vs_Sub(outputs, invData, bestRates, subsectionData, false, params, targets);

%% OPTIONAL TOOLS / EXPORTS
% plotLocalMFD(bestRates, invData, targets, subsectionData, params)
% exportInputToOpenQuake(bestRates, outputs, invData, All_Rupture_Candidates, faultResults, params)
% save_all_ruptures_SD_Batched(All_Final_Ruptures, subsectionData)

diary off;
fprintf('📝 Command Window log saved to: %s\n', logFileName);

%% 9. AUTOMATICALLY SAVE ALL OPEN FIGURES (PNG & SVG)
outputFolder = fullfile('Output', 'Plot');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('📁 Created directory: %s\n', outputFolder);
end

% Retrieve handles of all open figures
openFigures = findobj('Type', 'figure');

if ~isempty(openFigures)
    fprintf('💾 Saving %d open figures to %s...\n', length(openFigures), outputFolder);
    
    for f = 1:length(openFigures)
        figHandle = openFigures(f);
        
        % Get figure name or generate a default one
        figName = get(figHandle, 'Name');
        if isempty(figName)
            figName = sprintf('Figure_%d', figHandle.Number);
        else
            % Sanitize string for valid filesystem characters
            figName = strrep(figName, ' ', '_');
            figName = regexprep(figName, '[^\w\.-]', '');
        end
        
        % Save PNG (Raster, 300 DPI)
        pngPath = fullfile(outputFolder, [figName, '.png']);
        exportgraphics(figHandle, pngPath, 'Resolution', 300);
        
        % Save SVG (Vector Graphics for high-quality publication)
        svgPath = fullfile(outputFolder, [figName, '.svg']);
        exportgraphics(figHandle, svgPath, 'ContentType', 'vector');
    end
    
    fprintf('✅ All figures successfully saved in PNG and SVG formats.\n');
else
    fprintf('ℹ️ No open figures found to save.\n');
end
