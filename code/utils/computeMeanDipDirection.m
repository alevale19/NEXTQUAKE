function meanDipDir = computeMeanDipDirection(coordsWGS)
% Computes the mean dip direction from a WGS-84 polyline.
% INPUT:
%   coordsWGS : [N×2] double matrix containing [lon, lat] in degrees.
% OUTPUT:
%   meanDipDir : scalar representing the mean dip direction (0° to 360°).

    % 1. Simplify the polyline geometry using Douglas-Peucker algorithm.
    reduced = reducepoly(coordsWGS, 0.01); 
    
    % Extract longitudes and latitudes from the simplified coordinates
    lons = reduced(:,1); % First column contains Longitude (X)
    lats = reduced(:,2); % Second column contains Latitude (Y)
    
    % 2. Vectorized azimuth calculation between consecutive points.
    azis = azimuth(lats(1:end-1), lons(1:end-1), lats(2:end), lons(2:end));
    
    % 3. Calculate the mean azimuth using circular statistics (directional mean).
    meanAzi = atan2d(mean(sind(azis)), mean(cosd(azis)));
    
    % 4. Apply Right-Hand Rule (RHR): Dip Direction is 90 degrees clockwise from strike/azimuth.
    % Ensure the final output is wrapped within the 0 to 360 degrees range.
    meanDipDir = mod(meanAzi + 90, 360);
end