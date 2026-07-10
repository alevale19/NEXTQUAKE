function Global_Subsections_Data = computeSubsectionCentroids(faultResults, timeSettings, params)
% COMPUTESUBSECTIONCENTROIDS - Calculates 3D centroids, metadata, and inherits
% time-dependent parameters (Telaps, Thist, CV) for all subsections based on 
% a normalized along-strike position mapping with 3-level fallback logic.
%
% Inputs:
%   faultResults : Structure array from processFaultDatabase
%   timeSettings : Table/Struct with columns: FaultName, StartFrac, EndFrac, Telaps, Thist, CV

    % 1. Pre-calculate total number of subsections to avoid O(n^2) dynamic growth
    total_subs = sum([faultResults.m] .* [faultResults.n]);
    
    % 2. Pre-allocate arrays for maximum performance
    Global_ID = (1:total_subs)';
    Fault_ID  = zeros(total_subs, 1);
    FaultName = cell(total_subs, 1);
    S_idx     = zeros(total_subs, 1);
    D_idx     = zeros(total_subs, 1);
    X_km      = zeros(total_subs, 1); Y_km  = zeros(total_subs, 1); Z_km  = zeros(total_subs, 1);
    X1_km     = zeros(total_subs, 1); Y1_km = zeros(total_subs, 1); Z1_km = zeros(total_subs, 1);
    X2_km     = zeros(total_subs, 1); Y2_km = zeros(total_subs, 1); Z2_km = zeros(total_subs, 1);
    
    % Added Native Geographic Outputs for OpenQuake seamless reading
    Lon1      = zeros(total_subs, 1); Lat1  = zeros(total_subs, 1);
    Lon2      = zeros(total_subs, 1); Lat2  = zeros(total_subs, 1);
    
    Azimuth   = zeros(total_subs, 1);
    Length_km = zeros(total_subs, 1); 
    Width_km  = zeros(total_subs, 1); 
    Area_km2  = zeros(total_subs, 1);
    Coupling  = zeros(total_subs, 1);
    
    % Time Dependent Columns inherited from user input file
    Telaps    = zeros(total_subs, 1);
    Thist     = zeros(total_subs, 1);
    CV        = zeros(total_subs, 1);
    
    curr_idx = 1;
    
    % Check if time settings are provided, otherwise initialize empty
    if nargin < 2 || isempty(timeSettings)
        timeSettings = table();
    end
    
    for nf = 1:length(faultResults)
        f = faultResults(nf);
        
        m_strike = f.m;
        n_dip    = f.n;
        dip_val  = f.Dip;
        z_top_fault_km = f.upperST;
        Delta_W_km = f.Width / n_dip;
        
        lat_mesh = f.lat_mesh;
        lon_mesh = f.lon_mesh;
        
        [E_mesh_all, N_mesh_all, ~] = deg2utm(lat_mesh(:), lon_mesh(:));
        E_mesh = reshape(E_mesh_all, size(lat_mesh));
        N_mesh = reshape(N_mesh_all, size(lon_mesh));
        
        % Filter user time settings specifically for this fault
        faultTimeRows = table();
        if ~isempty(timeSettings) && ismember('FaultName', timeSettings.Properties.VariableNames)
            matchIdx = strcmpi(timeSettings.FaultName, f.Name);
            faultTimeRows = timeSettings(matchIdx, :);
        end
        
        % Nested loop to navigate the grid (s: strike, d: dip)
        for s = 1:m_strike
            
            % --- CRUCIAL: NORMALIZED POSITION ALONG STRIKE (0.0 to 1.0) ---
            current_frac = (s - 0.5) / m_strike;
            
            % Initialize fallbacks
            user_Telaps = NaN;
            user_Thist  = NaN;
            user_CV     = NaN;
            
            % Scan the user file rows to see if current_frac falls within a defined segment
            if ~isempty(faultTimeRows)
                for r = 1:height(faultTimeRows)
                    if current_frac >= faultTimeRows.StartFrac(r) && current_frac <= faultTimeRows.EndFrac(r)
                        user_Telaps = params.currentYear - faultTimeRows.LastEventYr(r);
                        user_Thist  = faultTimeRows.Thist(r);
                        user_CV     = faultTimeRows.CV(r);
                        break; % Take the first valid match
                    end
                end
            end
            
            % --- UCERF3 FALLBACK FOR CV (Aperiodicity) IF NaN ---
            if isnan(user_CV)
                % Estimate expected Magnitude based on Total Fault Area to assign CV class
                total_fArea = f.Width * f.Length;
                expected_Mw = seismicScalingEngine(total_fArea, 'Area2Mw', params.faultMechanism);
                
                if expected_Mw >= 6.3 && expected_Mw <= 6.7
                    assigned_CV = 0.5;
                elseif expected_Mw > 6.7 && expected_Mw <= 7.2
                    assigned_CV = 0.4;
                elseif expected_Mw > 7.2 && expected_Mw <= 7.7
                    assigned_CV = 0.3;
                elseif expected_Mw > 7.7
                    assigned_CV = 0.2;
                else
                    assigned_CV = 0.5; % Default standard
                end
            else
                assigned_CV = user_CV;
            end
            
            for d = 1:n_dip
                
                % --- A) CENTROID COORDINATES ---
                X_avg_m = (E_mesh(s,d) + E_mesh(s+1,d) + E_mesh(s,d+1) + E_mesh(s+1,d+1)) / 4;
                Y_avg_m = (N_mesh(s,d) + N_mesh(s+1,d) + N_mesh(s,d+1) + N_mesh(s+1,d+1)) / 4;
                Z_centroid_km = z_top_fault_km + ((d - 0.5) * Delta_W_km) * sind(dip_val);
                
                % --- B) TOP EDGE COORDINATES (Corners 1 and 2) ---
                X1_m = E_mesh(s, d);   Y1_m = N_mesh(s, d);
                X2_m = E_mesh(s+1, d); Y2_m = N_mesh(s+1, d);
                Z_edge_km = z_top_fault_km + ((d - 1) * Delta_W_km) * sind(dip_val);
                
                % --- C) NATIVE GEOGRAPHIC EXTRACTION ---
                Lon1(curr_idx) = lon_mesh(s, d);
                Lat1(curr_idx) = lat_mesh(s, d);
                Lon2(curr_idx) = lon_mesh(s+1, d);
                Lat2(curr_idx) = lat_mesh(s+1, d);
                
                % --- D) AZIMUTH CALCULATION ---
                az = atan2d(X2_m - X1_m, Y2_m - Y1_m);
                if az < 0, az = az + 360; end
                
                % --- E) DATA ASSIGNMENT ---
                Fault_ID(curr_idx)  = nf;
                FaultName{curr_idx} = f.Name;
                S_idx(curr_idx)     = s;
                D_idx(curr_idx)     = d;
                
                X_km(curr_idx)      = X_avg_m / 1000;
                Y_km(curr_idx)      = Y_avg_m / 1000;
                Z_km(curr_idx)      = Z_centroid_km;
                
                X1_km(curr_idx)     = X1_m / 1000;
                Y1_km(curr_idx)     = Y1_m / 1000;
                Z1_km(curr_idx)     = Z_edge_km;
                
                X2_km(curr_idx)     = X2_m / 1000;
                Y2_km(curr_idx)     = Y2_m / 1000;
                Z2_km(curr_idx)     = Z_edge_km;
                
                Azimuth(curr_idx)   = az;
                
                cell_length_m  = hypot(X2_m - X1_m, Y2_m - Y1_m);
                cell_length_km = cell_length_m / 1000;
                cell_width_km  = Delta_W_km;
                cell_area_km2  = cell_length_km * cell_width_km;
                
                Length_km(curr_idx) = cell_length_km;
                Width_km(curr_idx)  = cell_width_km;
                Area_km2(curr_idx)  = cell_area_km2;
                Coupling(curr_idx)  = f.CouplingCoefficient;
                
                % Assign Time Parameters (Inherited along strike segment mapping)
                Telaps(curr_idx)    = user_Telaps;
                Thist(curr_idx)     = user_Thist;
                CV(curr_idx)        = assigned_CV;
                
                curr_idx = curr_idx + 1;
            end
        end
    end
    
    % 3. Create the final unified table
    Global_Subsections_Data = table(Global_ID, Fault_ID, FaultName, S_idx, D_idx, ...
        X_km, Y_km, Z_km, X1_km, Y1_km, Z1_km, X2_km, Y2_km, Z2_km, ...
        Lon1, Lat1, Lon2, Lat2, Azimuth, Length_km, Width_km, Area_km2, Coupling, ...
        Telaps, Thist, CV); 
    
    fprintf('✅ Step 1 finished: %d subsections generated for %d faults (Time settings mapped successfully).\n', ...
            total_subs, length(faultResults));
end

