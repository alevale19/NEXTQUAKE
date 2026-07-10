function All_Elementary_Ruptures = generateElementaryRuptures(faultResults, Global_Subsections_Data, params)
% GENERATEELEMENTARYRUPTURES - Generates all possible rectangular ruptures for each fault.
%
% Optimized to avoid table growth in loops and handle O(n^2) complexity within faults.

    % Load user parameters
    AR_MAX = params.maxAspectRatio;    
    MIN_SUBSECTIONS = params.minSubs;   
    
    % Pre-allocate a cell array to store results (much faster than table growth)
    % We estimate a large enough size and then trim it 
    max_est = 0;
    for F = faultResults
        max_est = max_est + (F.m*(F.m+1)/2) * (F.n*(F.n+1)/2);
    end
    
    temp_results = cell(max_est, 10);
    rupt_count = 0;
    discarded_count = 0;

    faultIDs = unique(Global_Subsections_Data.Fault_ID);
    
    fprintf('🧱 Generating elementary ruptures for %d faults...\n', length(faultIDs));

    for f = 1:length(faultIDs)
        fid = faultIDs(f);
        
        % Filter data for current fault
        fData = Global_Subsections_Data(Global_Subsections_Data.Fault_ID == fid, :);
        m_strike = max(fData.S_idx);
        n_dip = max(fData.D_idx);
        
        % --- ROBUST CELL DIMENSION CALCULATION ---
        % Delta_L: Average distance between adjacent strike nodes
        if m_strike > 1
            Delta_L = faultResults(fid).Length / m_strike;
        else
            % If m_strike is 1, the subsection length is the full fault length
            Delta_L = faultResults(fid).Length;
        end
        
        % Delta_W: Average width along dip
        Delta_W = faultResults(fid).Width / max(1, n_dip);

        % --- QUAD-LOOP (Optimized) ---
        for s1 = 1:m_strike
            for s2 = s1:m_strike
                Ns = (s2 - s1 + 1);
                Length = Ns * Delta_L;
                
                for d1 = 1:n_dip
                    for d2 = d1:n_dip
                        Nd = (d2 - d1 + 1);
                        Width = Nd * Delta_W;
                        
                        n_subs = Ns * Nd;
                        
                        % Filter 1: Minimum size
                        if n_subs < MIN_SUBSECTIONS, continue; end
                        
                        % Filter 2: Aspect Ratio (L/W or W/L)
                        AR = max(Length/Width, Width/Length);
                        if AR > AR_MAX
                            discarded_count = discarded_count + 1; 
                            continue; 
                        end
                        
                        % Extract Global IDs for this rectangle
                        % Optimized logical indexing
                        mask = (fData.S_idx >= s1 & fData.S_idx <= s2 & ...
                                fData.D_idx >= d1 & fData.D_idx <= d2);
                        gIDs = fData.Global_ID(mask);

                        % Calculate Geometric Properties
                        Area = sum(fData.Area_km2(mask));
                        
                        Mw = seismicScalingEngine(Area, 'Area2Mw', params.faultMechanism);

                        % Calculate Mean Coordinates (Center of Rupture)
                        meanX = mean(fData.X_km(mask));
                        meanY = mean(fData.Y_km(mask));
                        meanZ = mean(fData.Z_km(mask));

                        % Store in temporary cell
                        rupt_count = rupt_count + 1;
                        temp_results(rupt_count, :) = {rupt_count, fid, s1, s2, d1, d2, ...
                                                       Area, Mw, gIDs, [meanX, meanY, meanZ]};
                    end
                end
            end
        end
    end

    % Convert to Table only at the end (Efficient)
    All_Elementary_Ruptures = cell2table(temp_results(1:rupt_count, :), ...
        'VariableNames', {'Rupture_ID', 'Fault_ID', 'S1', 'S2', 'D1', 'D2', ...
                          'Area', 'Mw', 'Global_IDs', 'Center_XYZ'});

    fprintf('✅ Generated %d elementary ruptures. (Discarded %d for Aspect Ratio)\n', rupt_count, discarded_count);
end