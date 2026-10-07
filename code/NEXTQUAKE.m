%% NEXTQUAKE - Main Processing Script
clc; clear; close all;
set(groot, 'DefaultFigureWindowStyle', 'docked');
addpath('./utils');

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
    
    % Load Constraints (Slip Rates, Paleo, MFD from files...)
    constraints = loadInversionConstraints(subsectionData, params);

    invData = prepareInversionData(All_Rupture_Candidates, subsectionData, constraints, params);
    
    % --- E. Save for Future Use ---
    save(invDataPath, 'invData', 'subsectionData', 'faultResults', 'params', 'A_S2S', '-v7.3');
    fprintf('💾 Pre-processed data saved to: %s\n', invDataPath);
end

%% 4. INVERSION CONSTRAINTS & TARGETS
% Setup Inversion Targets
% We merge the geometric data (invData) with physical constraints

targets = setupInversionTargets(invData, subsectionData, All_Rupture_Candidates, constraints, params);

%% 5. INITIALIZE RATES & TEST ENERGY
% Initialize all rupture rates to zero
% rates is a vector [numRuptures x 1]
rates = zeros(invData.numRuptures, 1);

[energy, components] = calculateEnergy(rates, invData, targets, params);

fprintf('Initial Energy: %.2e (SR: %.2e, MFD: %.2e, Paleo: %.2e, FaultMFD: %.2e, SegConstraint: %.2e)\n', ...
        energy, components.SR, components.MFD, components.Paleo, components.LocalMFD, components.Seg);

%% 6. LAUNCH SIMULATED ANNEALING
% --- BRIDGE: Mapping subsections and calculating total fault area ---
for i = 1:length(faultResults)
    % Find rows in subsectionData that match the current fault name
    matchedRows = strcmp(subsectionData.FaultName, faultResults(i).Name);
    
    % Store global IDs for the inversion logic
    faultResults(i).GlobalSubIndices = subsectionData.Global_ID(matchedRows);
    
    % Sum the area of all subsections to get the total fault area
    faultResults(i).TotalArea_km2 = sum(subsectionData.Area_km2(matchedRows));
end

[bestRates, energyHistory] = runSimulatedAnnealing(invData, targets, params, faultResults);

%% 7. QUICK RESULTS CHECK

% Fig.4 SA Convergence Plot, Fig.5 Event Rate match for each subsection with
% paleseismological rates (if any).

% Fig. 4
figure('Name', 'Convergence History');
semilogy(energyHistory, 'LineWidth', 1.5);
grid on;
xlabel('Iteration'); ylabel('Total Energy');
title('SA Convergence');

% Fig. 5
% --- Event Rate Match (Paleoseismic Comparison) ---
% 1. Calculate the model event rates for each subsection
finalEventRates = sum(invData.Gsr' .* bestRates.mean', 2);

figure('Name', 'Event Rate Match');
hold on;

% 2. Plot Paleo Sites as the primary reference
if isfield(invData, 'paleoTargetIdx') && ~isempty(invData.paleoTargetIdx)
    % Extract targets for clarity
    pIdx = invData.paleoTargetIdx;
    pRate = targets.paleoRate;
    pStd = targets.paleoStd;
    
    % Plot Error Bars for Paleo Sites (1-sigma)
    errorbar(pIdx, pRate, pStd, 'go', 'MarkerFaceColor', 'g', ...
        'LineWidth', 1.2, 'DisplayName', 'Paleo Targets (1\sigma)');
    
    % Plot the model results only for the subsections with paleo data (as dots)
    plot(pIdx, finalEventRates(pIdx), 'rs', 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Model at Paleo Sites');
end

% 3. Plot the continuous model event rate for all subsections
plot(finalEventRates, 'r-', 'LineWidth', 1, 'DisplayName', 'Model (mean)');

% Formatting
legend('show', 'Location', 'northeast');
grid on;
title('Event Rate Inversion Results & Paleoseismic Match');
ylabel('Event Rate (events/yr)');
xlabel('Subsection Index');

% Set Y-axis to Log scale if rates vary by orders of magnitude
set(gca, 'YScale', 'log');
%% 8. VISUALIZATION & QUALITY CHECK
fprintf('\n📊 Launching Visualizations...\n');

plotFaultSystem(faultResults, params); % Fig.6 2D Map Geometry & Subsections; Fig.7 Fault Clusters.
plotCentroids3D(faultResults, subsectionData); % Fig. 1 Plots 3D Centroids and Surface Fault Traces.
plotS2SConnections(faultResults, subsectionData, A_S2S); % Fig. 1 Subsection-to-Subsection (S2S) propagation bridges.
plotMFDComparison(bestRates, invData, targets, params); % Fig 1. & Fig.2 Plots Model MFD vs Target MFD with 5th and 95th Percentile Bounds (incremental and cumulative, respectively).

%% TOOLS

% === Calculates, formats, and saves standalone MFD plots for each individual fault === %
%plotLocalMFD(bestRates, invData, targets, subsectionData, params)

% === Compare the tectonic moment budget vs. the model release === %
checkMomentBalance(bestRates, invData, subsectionData, targets, params) % Fig.1 moment budget comparison

% === Spatial hazard visualizations optimized for 3D meshes === %
% 1. Fault Slip Rates (Constraint), 2. Mean Return Periods (Parent Fault Level), 3. Poisson Probabilities (Parent Fault Level)
% 4. BPT Time-Dependent Probabilities (Parent Fault Level), 5. Probability Gain Map (3D Subsection Patch Level), 6. Participation Rate Map (3D Subsection Patch Level in Log10 Scale)
outputs = processAnnealingResults(bestRates, invData, faultResults, subsectionData, params);
plotInversionAnalysisMaps(outputs, faultResults, subsectionData, invData, bestRates, params)

% === Analyzes and visualizes how frequently parent faults rupture together in multi-fault inversion solutions === %
analyzeFaultCoRupturePairs(faultResults, subsectionData, invData, bestRates)
analyzeComplexMultiFaultClusters(faultResults, subsectionData, invData, bestRates)
reportMultiFaultContribution(faultResults, subsectionData, invData, bestRates)

% === Plot Participation Rate at the subsection level. True for BPT-scaled rates, False for Poisson rates === %
plotGlobalParticipationM_vs_Sub(outputs, invData, bestRates, subsectionData, false, params, targets)

% === Save OpenQuake Input File === %
% exportInputToOpenQuake(bestRates, outputs, invData, All_Rupture_Candidates, faultResults, params)

% === Save plots for each Rupture === %
% save_all_ruptures_SD_Batched(All_Final_Ruptures, subsectionData)
