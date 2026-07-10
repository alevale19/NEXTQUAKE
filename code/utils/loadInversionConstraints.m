function constraints = loadInversionConstraints(subsectionData, params)
% LOADINVERSIONCONSTRAINTS - Loads Slip Rates, Paleo, Local MFD, and Segmentation rules.
% Mapping is performed using Euclidean distance in km (X,Y) and fault name matching.

fprintf('📥 Loading Inversion Constraints...\n');
numSubs = height(subsectionData);

% Initialize main structures
constraints.minSR = nan(numSubs, 1);     % m/yr
constraints.maxSR = nan(numSubs, 1);     % m/yr
constraints.stdSR = nan(numSubs, 1);     % m/yr
constraints.paleoRate = nan(numSubs, 1);   % event/yr
constraints.paleoStd = nan(numSubs, 1);    % event/yr
constraints.localBValue = ones(numSubs, 1);
constraints.badRupturePairs = [];          % List of forbidden fault name pairs

%% 1. MANDATORY: Fault-based Slip Rate (Box-Constraint)
if exist(params.fileSR, 'file')
    srTable = readtable(params.fileSR);
    allFaultNames = unique(subsectionData.FaultName);
    
    for i = 1:numel(allFaultNames)
        fName = allFaultNames{i};
        idx_table = find(strcmpi(srTable.FaultName, fName), 1);
        
        if ~isempty(idx_table)
        sub_idxs = find(strcmpi(subsectionData.FaultName, fName));

            % Convert mm/yr to m/yr
            sMin = (srTable.SRmin(idx_table) / 1000);
            sMax = (srTable.SRmax(idx_table) / 1000);
            
            constraints.minSR(sub_idxs) = sMin;
            constraints.maxSR(sub_idxs) = sMax;
            constraints.stdSR(sub_idxs) = (sMax - sMin) / 2;
        end
    end
else
    error('❌ Mandatory Slip Rate file not found: %s', params.fileSR);
end

%% 2. OPTIONAL: Point-based Slip Rate (Local corrections)
if isfield(params, 'filePointSR') && exist(params.filePointSR, 'file')
    % Load point slip rate data using headers from the .txt file
    pSR = readtable(params.filePointSR); 
    
    % Initialize counter to track successfully mapped point slip rates
    mappedPointCount = 0;
    numPointSites = height(pSR);

    for i = 1:numPointSites
        % Project geographic coordinates to the local Cartesian system (km)
        [pX, pY] = projectCoordinates(pSR.Lat(i), pSR.Lon(i), params);
        
        % Calculate distance from the point site to all subsections
        dist = sqrt((subsectionData.X_km - pX).^2 + (subsectionData.Y_km - pY).^2);
        
        % Filter by Fault Name to ensure the point is assigned to the correct fault
        faultMatch = strcmpi(subsectionData.FaultName, pSR.FaultName{i});
        
        if any(faultMatch)
            % If the fault name exists in the model, ignore distances to other faults
            dist(~faultMatch) = inf;
            [~, targetIdx] = min(dist);
            
            % Assign localized slip rate bounds (converting mm/yr to m/yr)
            constraints.minSR(targetIdx) = pSR.SRmin(i) / 1000;
            constraints.maxSR(targetIdx) = pSR.SRmax(i) / 1000;
            
            % Update standard deviation based on the new local bounds
            constraints.stdSR(targetIdx) = (constraints.maxSR(targetIdx) - constraints.minSR(targetIdx)) / 2;
            
            mappedPointCount = mappedPointCount + 1;
            fprintf('📍 Point SR [%s] mapped to Subsection #%d (Name match found).\n', ...
                pSR.FaultName{i}, targetIdx);
        else
            % If no subsections in the model match the FaultName in the point SR file
            fprintf('❌ WARNING: Point SR site [%s] could not be mapped because the Fault Name was not found in the model.\n', ...
                pSR.FaultName{i});
        end
    end
    
    % Final check: warn the user if any point-based rates were missed
    if mappedPointCount < numPointSites
        fprintf('\n⚠️  INCOMPLETE POINT SR MAPPING: Only %d out of %d sites were assigned.\n', ...
            mappedPointCount, numPointSites);
        fprintf('   Ensure the FaultName in the point SR file matches the subsection geometry exactly.\n\n');
    else
        fprintf('✅ Success: All %d point slip rate sites mapped successfully.\n', numPointSites);
    end
end

%% 3. OPTIONAL: Paleoseismic Event Rates
if isfield(params, 'filePaleo') && exist(params.filePaleo, 'file')
    % Read the table using the headers already present in the .txt file
    pal = readtable(params.filePaleo); 
    
    % Initialize a counter to track successfully mapped sites
    mappedCount = 0;
    numPaleoSites = height(pal);

    for i = 1:numPaleoSites
        % Project geographic coordinates to the local Cartesian system (km)
        [pX, pY] = projectCoordinates(pal.Lat(i), pal.Lon(i), params);
        
        % Calculate distance from the paleo site to all subsections
        dist = sqrt((subsectionData.X_km - pX).^2 + (subsectionData.Y_km - pY).^2);
        
        % Filter by Fault Name: prioritize finding the closest subsection 
        % that belongs to the specific fault mentioned in the paleo file.
        faultMatch = strcmpi(subsectionData.FaultName, pal.FaultName{i});
        
        if any(faultMatch)
            % If the fault name exists in the model, ignore distances to other faults
            dist(~faultMatch) = inf;
            [~, targetIdx] = min(dist);
            
            % Assign the Mean Rate to the closest subsection of that fault
            constraints.paleoRate(targetIdx) = pal.MeanRate(i);
            
            % Calculate 1-sigma standard deviation from the 95% confidence interval
            % Assuming pal.MinRate is 2.5% and pal.MaxRate is 97.5%
            lowerBound = pal.MinRate(i);
            upperBound = pal.MaxRate(i);
            constraints.paleoStd(targetIdx) = (upperBound - lowerBound) / 3.92;
            
            mappedCount = mappedCount + 1;
            fprintf('🏺 Paleo site [%s] mapped to Subsection #%d (Name match found).\n', ...
                pal.FaultName{i}, targetIdx);
        else
            % If no subsections in the model match the FaultName in the paleo file
            fprintf('❌ WARNING: Paleo site [%s] could not be mapped because the Fault Name was not found in the model.\n', ...
                pal.FaultName{i});
        end
    end
    
    % Final check: warn the user if the number of mapped sites is less than the file content
    if mappedCount < numPaleoSites
        fprintf('\n⚠️  INCOMPLETE PALEO MAPPING: Only %d out of %d sites were assigned.\n', ...
            mappedCount, numPaleoSites);
        fprintf('   Check for typos in Fault Names within your paleo .txt file.\n\n');
    else
        fprintf('✅ Success: All %d paleo sites mapped successfully.\n', numPaleoSites);
    end
end

%% 4. OPTIONAL: Local MFD (b-value per fault)
if isfield(params, 'fileLocalMFD') && exist(params.fileLocalMFD, 'file')
    bTable = readtable(params.fileLocalMFD); % Format: FaultName, bValue
    
    bTable.FaultName = strtrim(bTable.FaultName);
    
    for i = 1:height(bTable)
        sub_idxs = find(strcmpi(strtrim(subsectionData.FaultName), bTable.FaultName{i}));
        
        if isempty(sub_idxs)
            fprintf('⚠️  WARNING (Local MFD): Fault "%s" not found in geometry database!\n', bTable.FaultName{i});
            continue;
        end
        
        constraints.localBValue(sub_idxs) = bTable.bValue(i);
    end
    fprintf('🎯 Custom b-value constraints loaded for %d faults.\n', height(bTable));
end

%% 5. OPTIONAL: Segmentation Constraints (Forbidden pairs)
if isfield(params, 'fileSeg') && exist(params.fileSeg, 'file')
    constraints.segmentationPairs = readtable(params.fileSeg); % Format: FaultA, FaultB
    fprintf('🚧 Segmentation rules loaded (%d pairs).\n', height(constraints.segmentationPairs));
end

%% 6. REGIONAL MFD
if isfield(params, 'fileMFD') && exist(params.fileMFD, 'file')
    constraints.mfdData = readtable(params.fileMFD);
else
    constraints.mfdData = [];
end

fprintf('✅ Constraints loading complete.\n');
end

function [X, Y] = projectCoordinates(lat, lon, params)
    % PROJECTCOORDINATES - Converts Lat/Lon to UTM km
    [utmE, utmN, ~] = deg2utm(lat, lon);
    
    % IMPORTANT: If your subsectionData.X_km/Y_km are relative to a 
    % specific local origin (e.g. min(X), min(Y)), subtract it here.
    % If they are absolute UTM km:
    X = utmE / 1000; 
    Y = utmN / 1000;
end