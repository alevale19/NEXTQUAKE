function outputs = processAnnealingResults(bestRates, invData, faultResults, subsectionData, params)
    fprintf('🧮 Processing Simulated Annealing results...\n');
    
    if nargin < 5, params = struct(); end
    if isfield(params, 'forTwind'), t = params.forTwind; else, t = 30; end
    if isfield(params, 'forecastMw'), mw_lim = params.forecastMw; else, mw_lim = 6; end
    
    numSubsections = height(subsectionData);
    numFaults = length(faultResults);
    numRuptures = length(bestRates.mean);
    logicalGsr = invData.Gsr > 0;
    
    % --- Step 1: Subsection Stationarity (Poisson) ---
    subEventRate = (bestRates.mean(:)' * logicalGsr)'; 
    subTmean = 1 ./ subEventRate;
    subTmean(subEventRate == 0) = NaN;
    
    subPoissonProb = 1 - exp(-t ./ subTmean);
    subPoissonProb(isnan(subTmean)) = 0;
    
    % --- Step 2: Core BPT Calculation directly on SUBSECTIONS ---
    subBPTProb = zeros(numSubsections, 1);
    
    id_warn = 'MATLAB:integral:MaxIntervalCountReached';
    warning('off', id_warn);
    warning('off', 'MATLAB:integral:NonFiniteValue');
    
    for s = 1:numSubsections
        if isnan(subTmean(s)) || subTmean(s) == 0, continue; end
        
        mu_s = subTmean(s); 
        cv_s = subsectionData.CV(s);
        if cv_s <= 0 || isnan(cv_s), cv_s = 0.5; end
        
        lambda_s = mu_s / (cv_s^2);
        
        % BPT Funcion
        pdfBPT = @(T) pdf('inversegaussian', T, mu_s, lambda_s);
        cdfBPT = @(T) cdf('inversegaussian', T, mu_s, lambda_s);
        
        % 3-LEVEL FALLBACK CASCADE LOGIC
        %1. Time elapsed known
        if ~isnan(subsectionData.Telaps(s))
            q1 = cdfBPT(subsectionData.Telaps(s));
            q2 = cdfBPT(subsectionData.Telaps(s) + t);
            if q1 < 0.9999
                subBPTProb(s) = (q2 - q1) / (1 - q1); 
            else
                subBPTProb(s) = 1.0; 
            end
        elseif ~isnan(subsectionData.Thist(s))
            % Time elapsed unknow - use historic open interval
            S_func = @(T) (1 - cdfBPT(T));
            denom = integral(S_func, subsectionData.Thist(s), inf);
            num = integral(S_func, subsectionData.Thist(s), subsectionData.Thist(s) + t);
            if denom > 0
                subBPTProb(s) = num / denom;
            else
                subBPTProb(s) = subPoissonProb(s);
            end
        else
            % Time unkown
            q_unk = integral(@(T) (1 - cdfBPT(T)), 0, t);                      
            subBPTProb(s) = q_unk / mu_s;
        end
        
        if subBPTProb(s) > 1, subBPTProb(s) = 1; end
        if subBPTProb(s) < 0 || isnan(subBPTProb(s)), subBPTProb(s) = 0; end
    end
    warning('on', id_warn);
    
    % --- Step 3: compute subsections' Gain ---
    subGain = ones(numSubsections, 1);
    validGainIdx = subPoissonProb > 0;
    subGain(validGainIdx) = subBPTProb(validGainIdx) ./ subPoissonProb(validGainIdx);
    
    % --- Step 4: compute rupture's Gain and modify rates ---
    rupGain = ones(numRuptures, 1);
    bestRates_BPT = bestRates.mean(:); 
    secAreas = subsectionData.Area_km2; 
    
    for r = 1:numRuptures
        if bestRates.mean(r) == 0, continue; end 
        
        subInvolved = find(logicalGsr(r, :)); 
        if isempty(subInvolved), continue; end
        
        areas = secAreas(subInvolved);
        gains = subGain(subInvolved);
        
        totalRupArea = sum(areas);
        if totalRupArea > 0
            rupGain(r) = sum(gains .* areas) / totalRupArea;
        end
        
        bestRates_BPT(r) = bestRates.mean(r) * rupGain(r);
    end
    
    % --- Step 5: Parent Fault level aggregation ---
    faultParticipationRate_Poisson = zeros(numFaults, 1);
    faultParticipationRate_BPT     = zeros(numFaults, 1);
    faultPoissonProb               = zeros(numFaults, 1);
    faultBPTProb                   = zeros(numFaults, 1);
    
    for nf = 1:numFaults
        subIdx = faultResults(nf).GlobalSubIndices;
        if isempty(subIdx), continue; end
        
        for r = 1:numRuptures
            if invData.Mw(r) >= mw_lim && any(logicalGsr(r, subIdx))
                faultParticipationRate_Poisson(nf) = faultParticipationRate_Poisson(nf) + bestRates.mean(r);
                faultParticipationRate_BPT(nf)     = faultParticipationRate_BPT(nf) + bestRates_BPT(r);
            end
        end
        
        faultPoissonProb(nf) = 1 - exp(-t * faultParticipationRate_Poisson(nf));
        faultBPTProb(nf)     = 1 - exp(-t * faultParticipationRate_BPT(nf));
    end
    
    % --- Step 6: Output ---
    outputs = struct();
    outputs.bestRates_BPT = bestRates_BPT; 
    outputs.RuptureGains = rupGain;
    
    outputs.Subsections = table((1:numSubsections)', subEventRate(:), subTmean(:), ...
        subPoissonProb(:), subBPTProb(:), subGain(:), ...
        'VariableNames', {'Subsection_ID', 'ParticipationRate', 'Tmean', 'PoissonProb', 'BPT_TimeDependentProb', 'Gain'});
        
    outputs.Faults = table({faultResults.Name}', faultParticipationRate_Poisson, faultParticipationRate_BPT, ...
        faultPoissonProb, faultBPTProb, ...
        'VariableNames', {'FaultName', 'ParticipationRate', 'BPTRate', 'PoissonProb', 'BPT_TimeDependentProb'});
end
