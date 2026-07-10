function invData = prepareInversionData(All_Rupture_Candidates, GSD, constraints, params)
% PREPAREINVERSIONDATA - Prepares matrices for SR, Paleo, Local MFD, and Segmentation.

fprintf('📊 Preparing inversion matrices for %d ruptures...\n', height(All_Rupture_Candidates));
numRuptures = height(All_Rupture_Candidates);
numSubsections = height(GSD);
mu = params.mu; 

% --- 1. Moment and Slip Calculation ---
Mo = 10.^(1.5 * All_Rupture_Candidates.Mw + 9.05); 
areas_m2 = All_Rupture_Candidates.Area * 1e6; 
Dr = Mo ./ (mu * areas_m2);

% --- 2. Gsr Matrix (Standard Rupture-to-Subsection mapping) ---
row_idx = cell(numRuptures, 1);
col_idx = cell(numRuptures, 1);
for i = 1:numRuptures
    subs = All_Rupture_Candidates.Global_IDs{i};
    row_idx{i} = repmat(i, numel(subs), 1);
    col_idx{i} = subs(:);
end
Gsr = sparse(vertcat(row_idx{:}), vertcat(col_idx{:}), 1, numRuptures, numSubsections);

% --- 3. Gpaleo Matrix (Paleoseismic Mapping) ---
% Extract columns of Gsr corresponding to subsections with paleo data
paleoActiveIdx = find(~isnan(constraints.paleoRate));

if ~isempty(paleoActiveIdx)
    % 1. Map the base geometric matrix first
    invData.Gpaleo = Gsr(:, paleoActiveIdx);
    invData.paleoTargetIdx = paleoActiveIdx; % Map back to constraints
    
    % 2. Apply the magnitude-dependent visibility scaling if enabled
    if isfield(params, 'applyProbVisible') && params.applyProbVisible == 1
        % Calculate the probability that each given rupture magnitude will breach the surface
        rupMag = All_Rupture_Candidates.Mw(:); 
        
% --- Logistic regression for surface rupture visibility based on kinematics ---
        % Vectors of coefficients corresponding to: [1: All/Otherwise/Sub, 2: N, 3: R, 4: SS]
        % Coefficients come from Pizza et al., 2023
        coeff_a = [-14.47; -13.50; -10.75; -28.56];
        coeff_b = [  2.177;   2.159;   1.427;   4.436];
        
        switch upper(params.faultMechanism)
            case 'N'
                kinIdx = 2; % Normal faulting relations
            case 'R'
                kinIdx = 3; % Reverse / Thrust relations
            case 'SS'
                kinIdx = 4; % Strike-Slip relations
            case 'SUB'
                kinIdx = 1; % Subduction Zones (Uses ALL coefficients)
            otherwise
                kinIdx = 1; % Crustal ensemble / All (Uses ALL coefficients)
        end
        
        % Extract the specific coefficients for the selected mechanism
        a_val = coeff_a(kinIdx);
        b_val = coeff_b(kinIdx);
        
        % Compute the logistic function over the rupture magnitudes vector
        f_val = a_val + b_val .* rupMag;
        ProbVisible = exp(f_val) ./ (1 + exp(f_val));
        % Multiply each rupture row of the Gpaleo matrix by its visibility probability
        invData.Gpaleo = bsxfun(@times, invData.Gpaleo, ProbVisible);
        
        fprintf('📜 Paleoseismic matrix Gpaleo updated with magnitude-dependent surface visibility probability.\n');
    else
        fprintf('📜 Paleoseismic matrix Gpaleo mapped using purely geometric criteria (Visibility Prob. Disabled).\n');
    end
else
    invData.Gpaleo = [];
end

% --- 4. Local MFD Indices (Grouping by Fault Name) ---
% For each fault in subsectionData, find all ruptures that involve it
allFaultNames = unique(GSD.FaultName);
invData.faultRuptureIndices = cell(numel(allFaultNames), 1);
for f = 1:numel(allFaultNames)
    % Find ruptures that contain at least one subsection of this fault
    % Optimization: use the pre-built Gsr matrix
    subIdxOfFault = find(strcmpi(GSD.FaultName, allFaultNames{f}));
    invData.faultRuptureIndices{f} = find(any(Gsr(:, subIdxOfFault), 2));
end

% --- 5. Segmentation: Identify "Bad" Ruptures ---
invData.badRupIdx = [];
if isfield(constraints, 'segmentationPairs') && ~isempty(constraints.segmentationPairs)
    isBad = false(numRuptures, 1);
    
    % Get the fault name for every subsection once for speed
    subFaultNames = GSD.FaultName;
    
    for p = 1:height(constraints.segmentationPairs)
        f1 = strtrim(constraints.segmentationPairs.FaultA{p});
        f2 = strtrim(constraints.segmentationPairs.FaultB{p});
        
        % Find global IDs of all subsections belonging to Fault 1 and Fault 2
        subsF1 = find(strcmpi(subFaultNames, f1));
        subsF2 = find(strcmpi(subFaultNames, f2));
        
        if isempty(subsF1) || isempty(subsF2)
            continue; % One of the faults in the pair doesn't exist in the geometry
        end
        
        % A rupture is "bad" if it has at least one subsection from F1 
        % AND at least one subsection from F2.
        % We check the Gsr matrix: any(Gsr(:, subsF1), 2) is a logical array 
        % showing which ruptures contain Fault 1.
        hasF1 = any(Gsr(:, subsF1), 2);
        hasF2 = any(Gsr(:, subsF2), 2);
        
        isBad = isBad | (hasF1 & hasF2);
    end
    
    invData.badRupIdx = find(isBad);
    fprintf('🚧 Identified %d ruptures violating segmentation rules.\n', numel(invData.badRupIdx));
end

% --- 6. Final Structure ---
invData.numRuptures = numRuptures;
invData.numSubsections = numSubsections;
invData.Gsr = Gsr;
invData.Dr = Dr;
invData.Mw = All_Rupture_Candidates.Mw;
invData.Mo = Mo;
invData.subLengths = GSD.Length_km;

numPaleo = 0;
    if isfield(invData, 'Gpaleo') && ~isempty(invData.Gpaleo)
        numPaleo = size(invData.Gpaleo, 2);
    end

    fprintf('✅ System ready. Gsr: [%d x %d], Gpaleo sites: %d\n', ...
        numRuptures, numSubsections, numPaleo);
end