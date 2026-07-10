function [cUTM, cWGS] = enforceRightHandRule(cUTM, cWGS, strike)

    az_trace = azimuth(cWGS(1,2), cWGS(1,1), cWGS(end,2), cWGS(end,1));
    
    d1 = abs(wrapTo180(az_trace - strike));
    d2 = abs(wrapTo180(az_trace - (strike + 180)));
    
    if d2 < d1
        cUTM = flipud(cUTM);
        cWGS = flipud(cWGS);
    end
end