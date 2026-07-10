function [coordsUTM, coordsWGS, utmZone] = prepareCoordinates(faultName, system)
% PREPARECOORDINATES  Load Lon/Lat txt and convert to UTM.
%   [utmXY,wgs,zone] = prepareCoordinates(faultName) reads the file
%   'MainFaults_coords/<faultName>.txt' assumed to contain two numeric
%   columns [lon  lat] in degrees (WGS-84). It returns:
%      coordsUTM – [N×2] Easting, Northing (m)
%      coordsWGS – [N×2] [lat lon] (deg)
%      utmZone   – char, e.g. '33 T'
%
%   The function warns and returns empty arrays if the file is missing.
%
%   Limitations:
%     • Assumes the whole fault lies in a single UTM zone;
%     • Not valid poleward of 84° N / 80° S (UTM not defined);
%     • Input file must be numeric Lon Lat without header.
    
    fullPath = fullfile('MainFaults_coords', [faultName, '.txt']);
    
    if ~isfile(fullPath)
        warning('Coordinates file missing for %s: %s', faultName, fullPath);
        coordsUTM = []; coordsWGS = []; utmZone = ''; 
        return;
    end
    
    % Load raw data (Expected: Longitude, Latitude)
    data = load(fullPath); 
    
    % Even if 'system' is now always WGS84, we keep the logic to 
    % perform the conversion to UTM which is essential for geometry
    if strcmpi(system, 'WGS84')
        coordsWGS = data;
        
        % CONVERSION: Grids to Meters
        % deg2utm returns Easting, Northing, and the UTM Zone (e.g., '33 T')
        [E, N, utmZoneRaw] = deg2utm(coordsWGS(:,2), coordsWGS(:,1));
        coordsUTM = [E, N];
        
        % Check if utmZoneRaw is a matrix (multiple identical rows)
        % and force it to a single row string
        if size(utmZoneRaw, 1) > 1
            utmZone = utmZoneRaw(1, :);
        else
            utmZone = utmZoneRaw;
        end
        
    else
        % Fallback case (if you ever decide to use UTM files directly)
        coordsUTM = data;
        coordsWGS = data; % This would need a utm2deg conversion to be accurate
        utmZone = ' 33T'; % Default placeholder
    end
    
    % Ensure utmZone is always 4 characters long for utm2deg compatibility
    utmZone = char(utmZone);
    utmZone = strtrim(utmZone);
    if length(utmZone) == 3
        utmZone = [' ' utmZone];
    end
end