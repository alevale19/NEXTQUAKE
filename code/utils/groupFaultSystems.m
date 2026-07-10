function faultResults = groupFaultSystems(faultResults, params)
% GROUPFAULTSYSTEMS - Clusters faults into systems based on geometric and kinematic criteria.
% Uses graph theory to find connected components efficiently.

    nFault = length(faultResults);
    
    % Retrieve control parameters from the params structure (integrated tolerances)
    maxGap      = params.maxGap;
    maxAzDiff   = params.maxAzDiff;
    maxRakeDiff = params.maxRakeDiff;
    maxDipDiff  = params.maxDipDiff;
    
    % Initialize the logical connectivity adjacency matrix
    canConnect = false(nFault);
    
    % Extract vectorized attributes to avoid structural array overhead within the loop
    % This significantly boosts the execution speed
    allAz    = [faultResults.Azimuth];
    allRake  = [faultResults.Rake];
    allDip   = [faultResults.Dip];
    
    % Prepare trace segment endpoints (Start-End)
    % Extract existing UTM coordinates and convert them to kilometers
    starts = zeros(nFault, 2);
    ends   = zeros(nFault, 2);
    for i = 1:nFault
        coords = faultResults(i).coordsUTM / 1000; % Conversion to km
        starts(i, :) = coords(1, :);
        ends(i, :)   = coords(end, :);
    end
    
    fprintf('🔍 Clustering %d faults into systems...\n', nFault);
    
    % --- OPTIMIZED CONNECTIVITY LOOP ---
    for i = 1:nFault
        for j = i+1:nFault % Compare each pair exactly once
            
            % 1. Minimum Distance Between Segments (Structural Gap Check)
            d = distBetweenSegments(starts(i,:), ends(i,:), starts(j,:), ends(j,:));
            if d > maxGap, continue; end % Early pruning if segments are too far apart
            
            % 2. Azimuth Difference (Handling 180-degree boundary circularity)
            dAz = abs(allAz(i) - allAz(j));
            dAz = mod(dAz + 90, 180) - 90;
            if abs(dAz) > maxAzDiff, continue; end
            
            % 3. Rake Difference (Kinematic Compatibility Check)
            dRake = abs(allRake(i) - allRake(j));
            if dRake > 180, dRake = 360 - dRake; end
            if dRake > maxRakeDiff, continue; end
            
            % 4. Dip Angle Difference (Geometric Coherency Check)
            dDip = abs(allDip(i) - allDip(j));
            if dDip > maxDipDiff, continue; end
            
            % If all structural filters pass, mark them as connected
            canConnect(i,j) = true;
            canConnect(j,i) = true;
        end
    end
    
    % --- COMPONENT COMPUTATION (Graph Theory) ---
    % Convert the logical adjacency matrix into a MATLAB graph object
    G = graph(canConnect);
    
    % Extract connected components (systemIDs(i) yields the group ID for fault i)
    systemIDs = conncomp(G);
    
    % --- ASSIGN RESULTS TO STRUCT ---
    for i = 1:nFault
        faultResults(i).SystemID = systemIDs(i);
    end
    
    nSystems = max(systemIDs);
    fprintf('✅ Clustering complete: Found %d independent fault systems.\n', nSystems);
end