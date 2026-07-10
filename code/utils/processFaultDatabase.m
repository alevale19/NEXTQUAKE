function results = processFaultDatabase(Faultdata, params, subSectLength)
% PROCESSFAULTDATABASE - Processes the fault database to build a 3D discretization mesh.
% This version handles custom upperST and lowerST bounds, projects the sismogenic grid
% starting from the upper seismogenic edge, and computes the real mean azimuth from nodes.
% It also includes a pre-check warning if the shortest fault restricts discretization.
%
% Inputs:
%   Faultdata      : Table containing fault parameters (FaultName, Dip, Strike, upperST, lowerST, Rake)
%   params         : Structure containing system configurations (e.g., coordSystem, forecastMw, faultMechanism)
%   subSectLength  : Desired subsection length along strike and dip width (in km)
%
% Outputs:
%   results        : Structure array containing full geometric properties and meshes

    n_Mainfault = size(Faultdata, 1);
    fprintf('You have %i Mainfaults in the DB\n', n_Mainfault);
   
    % =================================================================
    % PRE-DISCRETIZATION DIAGNOSTIC & WARNING
    % =================================================================
    % Check lengths of all faults in the database to catch sub-discretization issues
    all_lengths_km = zeros(n_Mainfault, 1);
    all_areas_km2  = zeros(n_Mainfault, 1);
    
    for nf = 1:n_Mainfault
        faultName = char(Faultdata.FaultName(nf));
        [coordsUTM, ~, ~] = prepareCoordinates(faultName, params.coordSystem);
        if ~isempty(coordsUTM)
            all_lengths_km(nf) = sum(hypot(diff(coordsUTM(:,1)), diff(coordsUTM(:,2)))) / 1000;
            faultWidth = (Faultdata.lowerST(nf) - Faultdata.upperST(nf)) / sind(Faultdata.Dip(nf));
            all_areas_km2(nf)  = all_lengths_km(nf) * faultWidth;
        end
    end
    
    % Find the critical (shortest) fault
    [minFaultLength, minIdx] = min(all_lengths_km(all_lengths_km > 0));
    criticalFaultName = char(Faultdata.FaultName(minIdx));
    
if minFaultLength < subSectLength
        % Calculate the physical magnitude bound of this critical fault
        criticalArea = all_areas_km2(minIdx);
        suggestedMw = seismicScalingEngine(criticalArea, 'Area2Mw', params.faultMechanism);
        
        % TRIGGER WARNING ONLY IF THE DATABASE TRULY UNDERESTIMATES THE TARGET MW
        if params.forecastMw > suggestedMw
            warning_msg = [ ...
                'The shortest fault trace in the database (%s: %.2f km) is shorter than the chosen subsection length (%.2f km).\n' ...
            'As a result, your target forecast Mw (%.2f) cannot be properly discretized on this fault.\n' ...
            'The actual magnitude floor will drop, meaning the forecast Mw could be underestimated.\n' ...
            '--> SUGGESTION: Update your forecast Mw to approximately %.2f (or lower) to match the database constraints.\n' ...
            '========================================================================================\n'];

            warning(warning_msg, criticalFaultName, minFaultLength, subSectLength, params.forecastMw, suggestedMw);
        end
end
    % =================================================================
    % MAIN PROCESSING LOOP
    % =================================================================
    results = struct(); 
    for nf = 1:n_Mainfault
        faultName = char(Faultdata.FaultName(nf));
        
        % 1) LOAD AND CONVERT COORDINATES
        [coordsUTM, coordsWGS, utmZone] = prepareCoordinates(faultName, params.coordSystem);
        if isempty(coordsUTM), continue; end
        
        % 2) RIGHT-HAND RULE (RHR) ENFORCEMENT
        strike_val = Faultdata.Strike(nf);
        [coordsUTM, coordsWGS] = enforceRightHandRule(coordsUTM, coordsWGS, strike_val);
        
        % 3) GEOMETRY PARAMETERS
        dip_deg = Faultdata.Dip(nf);
        seisThick = (Faultdata.lowerST(nf) - Faultdata.upperST(nf)); 
        rake = Faultdata.Rake(nf);
        faultWidth = seisThick / sind(dip_deg); 
        
        % 3.a) CREEP TREATMENT
        creep_fraction = Faultdata.CreepFraction(nf);  
        if creep_fraction <= 0.4
            aseismicity_factor = creep_fraction;
            coupling_coefficient = 1.0;
        else
            aseismicity_factor = 0.4;
            coupling_coefficient = 1.0 - ((creep_fraction - 0.4) / (1.0 - 0.4));
        end
        
        % 4) FAULT LENGTH & DISCRETIZATION SEGMENTS
        totalLength_km = all_lengths_km(nf);
        m_segments = max(floor(totalLength_km / subSectLength), 1); 
        n_segments = max(floor(faultWidth / subSectLength), 1);    
        
        % 5) DIP DIRECTION
        meanDipDir = computeMeanDipDirection(coordsWGS);
        direction_rad = deg2rad(meanDipDir);
        
        % 6) UTM ZONE SANITIZATION & DEPTH PROJECTIONS
        if iscell(utmZone)
            uz = utmZone{1}; 
        else
            uz = utmZone;
        end
        uz = char(uz); uz = uz(1, :); uz = strtrim(uz);
        if length(uz) == 3
            uz = [' ' uz]; 
        elseif length(uz) > 4
            uz = uz(1:4);
        end
        results(nf).utmZone = uz;
        
        dist_upper_m = (Faultdata.upperST(nf) * 1000) / tand(dip_deg);
        delta_E_up = dist_upper_m * sin(direction_rad);
        delta_N_up = dist_upper_m * cos(direction_rad);
        
        dist_lower_m = (Faultdata.lowerST(nf) * 1000) / tand(dip_deg); 
        delta_E_low = dist_lower_m * sin(direction_rad);
        delta_N_low = dist_lower_m * cos(direction_rad);
        
        E_new = coordsUTM(:,1) + delta_E_low;
        N_new = coordsUTM(:,2) + delta_N_low;
        E_new = E_new(:); N_new = N_new(:);
       
        n_trace = length(E_new);
        z_matrix_trace = repmat(uz, n_trace, 1);
        [lat_new, lon_new] = utm2deg(E_new, N_new, z_matrix_trace);
        
        lat = coordsWGS(:,2); lon = coordsWGS(:,1);
        box = [lat, lon; flipud(lat_new), flipud(lon_new); lat(1), lon(1)];
        
        % 7) SUBSECTION NODES ALONG STRIKE (Spline Interpolation)
        distances = [0; cumsum(hypot(diff(coordsUTM(:,1)), diff(coordsUTM(:,2))))];
        target_distances = linspace(0, distances(end), m_segments + 1);
        
        E_nodes_top = spline(distances, coordsUTM(:,1), target_distances)';
        N_nodes_top = spline(distances, coordsUTM(:,2), target_distances)';
        n_nodes = length(E_nodes_top);
        z_matrix_nodes = repmat(uz, n_nodes, 1);
        
        [lat_interp, lon_interp] = utm2deg(E_nodes_top, N_nodes_top, z_matrix_nodes);
        
        E_trans = E_nodes_top + delta_E_low;
        N_trans = N_nodes_top + delta_N_low;
        [lat_translated, lon_translated] = utm2deg(E_trans(:), N_trans(:), z_matrix_nodes);
        
        % COMPUTE MEAN REAL GEOMETRIC AZIMUTH FROM NODES
        seg_dX = diff(E_nodes_top);
        seg_dY = diff(N_nodes_top);
        seg_azimuths = atan2d(seg_dX, seg_dY);
        seg_azimuths(seg_azimuths < 0) = seg_azimuths(seg_azimuths < 0) + 360;
        
        mean_sin = mean(sin(deg2rad(seg_azimuths)));
        mean_cos = mean(cos(deg2rad(seg_azimuths)));
        calculated_mean_azimuth = atan2d(mean_sin, mean_cos);
        if calculated_mean_azimuth < 0
            calculated_mean_azimuth = calculated_mean_azimuth + 360;
        end
        
        % 8) MESH DISCRETIZATION ALONG WIDTH (Interpolating Upper & Lower Edges)
        lat_mesh = zeros(m_segments+1, n_segments+1);
        lon_mesh = zeros(m_segments+1, n_segments+1);
        depth_mesh = zeros(m_segments+1, n_segments+1); 
        
        for j = 0:n_segments
            frac = j / n_segments;           
            E_current = E_nodes_top + (delta_E_up + (delta_E_low - delta_E_up) * frac);
            N_current = N_nodes_top + (delta_N_up + (delta_N_low - delta_N_up) * frac);         
            current_depth_km = Faultdata.upperST(nf) + (seisThick * frac);
            depth_mesh(:, j+1) = current_depth_km; 
            
            n_pts = length(E_current);
            z_mat = repmat(uz, n_pts, 1);
            [lats, lons] = utm2deg(E_current(:), N_current(:), z_mat);
            lat_mesh(:, j+1) = lats;
            lon_mesh(:, j+1) = lons;
        end
        
        % 9) SAVE EVERYTHING TO RESULTS STRUCT
        if isprop(Faultdata, 'Fault_ID') || ismember('Fault_ID', Faultdata.Properties.VariableNames)
            results(nf).Fault_ID = Faultdata.Fault_ID(nf);
        else
            results(nf).Fault_ID = nf; 
        end
        
        results(nf).Name = faultName;
        results(nf).Length = totalLength_km;
        results(nf).Width = faultWidth;
        results(nf).lowerST = Faultdata.lowerST(nf);
        results(nf).upperST = Faultdata.upperST(nf);
        results(nf).Dip = dip_deg;
        results(nf).Rake = rake;
        results(nf).Azimuth = calculated_mean_azimuth; 
        results(nf).AseismicityFactor = aseismicity_factor;
        results(nf).CouplingCoefficient = coupling_coefficient;
        
        results(nf).m = m_segments;
        results(nf).n = n_segments;
        results(nf).coordsWGS = coordsWGS;
        results(nf).coordsWGS_new = [lon_new, lat_new];
        results(nf).coordsUTM = coordsUTM;
        results(nf).lat_interp = lat_interp;
        results(nf).lon_interp = lon_interp;
        results(nf).lat_translated = lat_translated;
        results(nf).lon_translated = lon_translated;
        results(nf).box = box;
        results(nf).lat_mesh = lat_mesh;
        results(nf).lon_mesh = lon_mesh;
        results(nf).depth_mesh = depth_mesh; 
        
        fprintf('Processed: %s (ID=%i, m=%i, n=%i) | Zone: %s\n', ...
            faultName, results(nf).Fault_ID, m_segments, n_segments, uz);
    end
end