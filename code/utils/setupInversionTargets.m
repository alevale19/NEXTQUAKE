function targets = setupInversionTargets(invData, subsectionData, All_Rupture_Candidates, constraints, params)
% SETUPINVERSIONTARGETS - Prepares target vectors, MFD bins, and weights for SA.
% This version includes Local MFD (b-values), Paleo rates, and Segmentation weights.

fprintf('🎯 Finalizing Inversion Targets...\n');

%% 1. SLIP RATE TARGETS (Mean, Min, Max)
targets.minSR    = constraints.minSR;
targets.maxSR    = constraints.maxSR;
targets.stdSR    = constraints.stdSR;

% Safety check: avoid division by zero or NaN
min_std = 0.00001; % 0.01 mm/yr
invalid_std = (targets.stdSR <= 0 | isnan(targets.stdSR));
if any(invalid_std)
    % Assign a 20% uncertainty if no std is provided
    targets.stdSR(invalid_std) = max((targets.minSR(invalid_std) + targets.maxSR(invalid_std))/2 * 0.2, min_std);
end

%% 2. PALEOSEISMIC TARGETS
% Initialize as empty so the inversion knows there are no paleo constraints
targets.paleoRate = [];
targets.paleoStd  = [];

% Check if paleoTargetIdx exists AND is not empty
if isfield(invData, 'paleoTargetIdx') && ~isempty(invData.paleoTargetIdx)
    
    % Only extract rates for subsections identified as paleo sites
    targets.paleoRate = constraints.paleoRate(invData.paleoTargetIdx);
    targets.paleoStd  = constraints.paleoStd(invData.paleoTargetIdx);
    
    % Safety check for Paleo Std (to avoid division by zero or NaN)
    invalid_p_std = (targets.paleoStd <= 0 | isnan(targets.paleoStd));
    if any(invalid_p_std)
        % Default to 20% uncertainty if std is missing or invalid
        targets.paleoStd(invalid_p_std) = targets.paleoRate(invalid_p_std) * 0.2;
    end
    
    fprintf('🎯 Paleo Targets set for %d sites.\n', numel(targets.paleoRate));
else
    fprintf('ℹ️ No Paleoseismic targets to set (optional file not provided).\n');
end

%% 3. LOCAL MFD TARGETS (b-values per subsection)
targets.localBValue = constraints.localBValue;

%% 4. REGIONAL MFD TARGET
if isfield(constraints, 'mfdData') && ~isempty(constraints.mfdData)
    % --- CASE A: CUSTOM MFD FILE EXISTS ---
    % Directly extract bins and rates from the input file without modifications
    fprintf('📊 Using Custom Regional MFD directly from input file columns.\n');
    
    targets.mfdBins = constraints.mfdData.magBin(:);   % Extract magnitude bins from file
    targets.mfdTarget = constraints.mfdData.rate(:); % Extract exact historical rates from file
    
else
    % --- CASE B: NO MFD FILE FOUND (Energy-Balanced Truncated G-R) ---
    fprintf('📊 No MFD file found. Computing Energy-Balanced Truncated G-R Target on Catalog Range...\n');
    
    % --- Step A: Setup Constants ---
    shear_modulus = params.mu;
    MoRateReduction = params.MoRateReduction;
    bValue = params.b_val;  % Target b-value from the GUI (default 1.0)
    mag_delta = params.mag_delta;        % Magnitude bin width
    
    % --- Step B: Calculate True Seismic Moment Rate from Subsections ---
    L_m = subsectionData.Length_km * 1000; % Convert subsection length to meters
    W_m = subsectionData.Width_km * 1000;  % Convert original full width to meters
    
    % Retrieve creep-related coefficients inherited from the fault database
    if ismember('Coupling', subsectionData.Properties.VariableNames)
        coupling = subsectionData.Coupling;
        
        % Reconstruct the Aseismicity Factor from the coupling logic:
        % Under NSHM rules, if coupling < 1.0, creep fraction exceeded 0.4,
        % meaning the maximum area reduction (aseismicFactor = 0.4) has been reached.
        aseismicFactor = zeros(size(coupling));
        aseismicFactor(coupling < 1.0) = 0.4; 
    else
        coupling = ones(size(L_m));
        aseismicFactor = zeros(size(L_m));
    end
    
    % Extract nominal slip rates from constraints (unaltered by coupling in Box-Constraints)
    secSlipRate_m_yr = (constraints.minSR + constraints.maxSR) / 2; 
    
    % Seismic Moment Rate Calculation incorporating USGS NSHM Creep Treatment:
    % - Seismogenic Area is reduced by: (1 - aseismicFactor)
    % - Seismic Slip Rate is reduced by: secSlipRate_m_yr * coupling
    % This correctly implements the fractional reduction of seismic moment energy.
    secMoRate = shear_modulus .* (L_m .* (W_m .* (1 - aseismicFactor))) .* (secSlipRate_m_yr .* coupling);
    totalMoRate = sum(secMoRate) * (1 - MoRateReduction);
    
    % --- Step C: Define Magnitude Range Exactly from Rupture Catalog ---
    minMag = min(All_Rupture_Candidates.Mw); 
    maxMag = max(All_Rupture_Candidates.Mw);
    
    targets.mfdBins = (minMag : mag_delta : maxMag)';
    
    % Hanks & Kanamori (1979) Mw to Moment (Nm) conversion parameters
    c_hk = 1.5; f_hk = 9.05;
    M_bins = 10.^(c_hk .* targets.mfdBins + f_hk);
    Mt = 10^(c_hk * minMag + f_hk);
    Mxp = 10^(c_hk * (maxMag + mag_delta) + f_hk); % Safety padding bin for Max Magnitude
    
    Beta = (2/3) * bValue;
    
   % --- Step D: Truncated Gutenberg-Richter Core Formulation ---
        % Force TruncGR to be a strict column vector
        TruncGR = (((Mt./M_bins).^Beta) - ((Mt./Mxp).^Beta)) ./ (1 - (Mt./Mxp).^Beta);
        Incremental = -diff([TruncGR;0]);
        Incremental_Morate = Incremental .* M_bins;
        Incremental_Morate_balanced = (Incremental_Morate .* totalMoRate) ./ sum(Incremental_Morate);
        
        % Derive the balanced incremental seismic recurrence rates
        cons_tassi_ind = (Incremental .* Incremental_Morate_balanced) ./ (Incremental_Morate);
  
        % --- Step E: Target Assignment ---
    % Map the generated recurrence rates directly to the target vector structure
    targets.mfdTarget = cons_tassi_ind(:);
    
    fprintf('  -> Catalog Range: Mw %.1f to %.1f\n', minMag, maxMag);
    fprintf('  -> Total Seismic Moment Rate (Creep Balanced): %.2e Nm/yr\n', totalMoRate);
    fprintf('  -> Truncated G-R Target built successfully.\n');
end

numBins = numel(targets.mfdBins);

%% 5. RUPTURE-TO-BIN MAPPING (Optimization)
numRuptures = invData.numRuptures;
targets.rupBinIdx = zeros(numRuptures, 1);
for i = 1:numRuptures
    [~, binIdx] = min(abs(targets.mfdBins - invData.Mw(i)));
    targets.rupBinIdx(i) = binIdx;
end

%% 6. INVERSION WEIGHTS (The "Knobs")
targets.weightSR       = params.weightSR;
targets.weightMFD      = params.weightMFD;
targets.weightL2       = params.weightL2;
targets.weightPaleo    = params.weightPaleo;
targets.weightLocalMFD = params.weightLocalMFD;
targets.weightSeg      = params.weightSeg;

fprintf('✅ Targets ready: %d Slip Targets, %d Paleo Sites, %d MFD Bins.\n', ...
        sum(targets.maxSR > 0), numel(targets.paleoRate), numBins);
end