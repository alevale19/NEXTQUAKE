function outputs = processAnnealingResults(bestRates, invData, faultResults, subsectionData, params)
    fprintf('🧮 Processing Simulated Annealing results (Rupture-Level Area-Weighted Gain for OpenQuake)...\n');
    
    if nargin < 5, params = struct(); end
    if isfield(params, 'forTwind'), t = params.forTwind; else, t = 30; end
    if isfield(params, 'forecastMw'), mw_lim = params.forecastMw; else, mw_lim = 6.5; end
    
    numSubsections = height(subsectionData);
    numFaults = length(faultResults);
    numRuptures = length(bestRates.mean);
    logicalGsr = invData.Gsr > 0; % Matrice binaria Ruoture x Sottosezioni
    
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
        funBPT = @(T) cdf('inversegaussian', T, mu_s, lambda_s);
        
        % 3-LEVEL FALLBACK CASCADE LOGIC
        if ~isnan(subsectionData.Telaps(s))
            q1 = funBPT(subsectionData.Telaps(s));
            q2 = funBPT(subsectionData.Telaps(s) + t);
            if q1 < 0.9999, subBPTProb(s) = (q2 - q1) / (1 - q1); else, subBPTProb(s) = 1.0; end
        elseif ~isnan(subsectionData.Thist(s))
            q1 = integral(funBPT, subsectionData.Thist(s), subsectionData.Thist(s) + t);
            funDenom = @(T) (1 - cdf('inversegaussian', T, mu_s, lambda_s));
            q2 = integral(funDenom, subsectionData.Thist(s), inf);      
            subBPTProb(s) = (t - q1) / q2;
        else
            q_unk = integral(funBPT, 0, t);                      
            subBPTProb(s) = (t - q_unk) / mu_s;
        end
        
        if subBPTProb(s) > 1, subBPTProb(s) = 1; end
        if subBPTProb(s) < 0 || isnan(subBPTProb(s)), subBPTProb(s) = 0; end
    end
    warning('on', id_warn);
    
    % --- Step 3: Calcolo del Gain di ogni Sottosezione ---
    subGain = ones(numSubsections, 1);
    for s = 1:numSubsections
        if subPoissonProb(s) > 0
            subGain(s) = subBPTProb(s) / subPoissonProb(s);
        end
    end
    
    % --- Step 4: Calcolo del Gain per ROTTURA e Modifica dei Tassi ---
    rupGain = ones(numRuptures, 1);
    bestRates_BPT = bestRates.mean(:); % Inizializziamo il vettore dei tassi modificati
    
    secAreas = subsectionData.Area_km2; % Vettore delle aree delle sottosezioni
    
    for r = 1:numRuptures
        if bestRates.mean(r) == 0, continue; end % Salta le rotture inattive
        
        % Trova gli indici delle sottosezioni coinvolte in questa rottura
        subInvolved = find(logicalGsr(r, :)); 
        if isempty(subInvolved), continue; end
        
        % Estrai aree e gain di queste sottosezioni
        areas = secAreas(subInvolved);
        gains = subGain(subInvolved);
        
        totalRupArea = sum(areas);
        if totalRupArea > 0
            % MEDIA PESATA SULL'AREA del Gain per la rottura r
            rupGain(r) = sum(gains .* areas) / totalRupArea;
        end
        
        % MODIFICA DEL TASSO POISSONIANO ORIGINALE
        bestRates_BPT(r) = bestRates.mean(r) * rupGain(r);
    end
    
    % --- Step 5: Aggregazione finale a livello di Parent Fault (per i report) ---
    faultParticipationRate_Poisson = zeros(numFaults, 1);
    faultParticipationRate_BPT     = zeros(numFaults, 1);
    faultPoissonProb               = zeros(numFaults, 1);
    faultBPTProb                   = zeros(numFaults, 1);
    
    globalLookup = struct('nf', cell(numSubsections, 1));
    current_global = 1;
    for nf = 1:numFaults
        for s_idx = 1:faultResults(nf).m
            for d_idx = 1:faultResults(nf).n
                globalLookup(current_global).nf = nf;
                current_global = current_global + 1;
            end
        end
    end
    
    for nf = 1:numFaults
        subIdx = find([globalLookup.nf] == nf);
        if isempty(subIdx), continue; end
        
        for r = 1:numRuptures
            if invData.Mw(r) >= mw_lim && any(logicalGsr(r, subIdx))
                faultParticipationRate_Poisson(nf) = faultParticipationRate_Poisson(nf) + bestRates.mean(r);
                faultParticipationRate_BPT(nf)     = faultParticipationRate_BPT(nf) + bestRates_BPT(r);
            end
        end
        faultPoissonProb(nf) = 1 - exp(-t * faultParticipationRate_Poisson(nf));
        faultBPTProb(nf)     = 1 - prod(1 - subBPTProb(subIdx));
    end
    
   % --- Step 6: Compilazione Output (Nomi variabili allineati al 100% con plotInversionAnalysisMaps) ---
    outputs = struct();
    outputs.bestRates_BPT = bestRates_BPT; % Passato all'esportatore XML per OpenQuake
    outputs.RuptureGains = rupGain;
    
    outputs.Subsections = table((1:numSubsections)', subEventRate(:), subTmean(:), ...
        subPoissonProb(:), subBPTProb(:), subGain(:), ...
        'VariableNames', {'Subsection_ID', 'ParticipationRate', 'Tmean', 'PoissonProb', 'BPT_TimeDependentProb', 'Gain'});
        
    outputs.Faults = table({faultResults.Name}', faultParticipationRate_Poisson, faultParticipationRate_BPT, ...
        faultPoissonProb, faultBPTProb, ...
        'VariableNames', {'FaultName', 'ParticipationRate', 'BPTRate', 'PoissonProb', 'BPT_TimeDependentProb'});
end
% function outputs = processAnnealingResults(bestRates, invData, faultResults, subsectionData, params)
%     fprintf('🧮 Processing Simulated Annealing results (Rigorous Subsection-Level BPT)...\n');
% 
%     if nargin < 5, params = struct(); end
%     if isfield(params, 'forTwind'), t = params.forTwind; else, t = 30; end
%     if isfield(params, 'forecastMw'), mw_lim = params.forecastMw; else, mw_lim = 6.5; end
% 
%     numSubsections = height(subsectionData);
%     numFaults = length(faultResults);
%     logicalGsr = invData.Gsr > 0; 
% 
%     % --- Step 1: Subsection Stationarity (Poisson) ---
%     subEventRate = (bestRates.mean(:)' * logicalGsr)'; 
%     subTmean = 1 ./ subEventRate;
%     subTmean(subEventRate == 0) = NaN;
% 
%     subPoissonProb = 1 - exp(-t ./ subTmean);
%     subPoissonProb(isnan(subTmean)) = 0;
% 
%     % --- Step 2: Core BPT Calculation directly on SUBSECTIONS ---
%     subBPTProb = zeros(numSubsections, 1);
% 
%     id_warn = 'MATLAB:integral:MaxIntervalCountReached';
%     warning('off', id_warn);
%     warning('off', 'MATLAB:integral:NonFiniteValue');
% 
%     for s = 1:numSubsections
%         % Se la sottosezione è inattiva (tasso zero), la sua probabilità è zero
%         if isnan(subTmean(s)) || subTmean(s) == 0
%             continue;
%         end
% 
%         mu_s = subTmean(s); % Tempo di ritorno reale della sottosezione (es. 312 anni)
%         cv_s = subsectionData.CV(s);
%         if cv_s <= 0 || isnan(cv_s), cv_s = 0.5; end
% 
%         % Parametro di forma lambda corretto per la Inverse Gaussian di MATLAB
%         lambda_s = mu_s / (cv_s^2);
%         funBPT = @(T) cdf('inversegaussian', T, mu_s, lambda_s);
% 
%         % 3-LEVEL FALLBACK CASCADE LOGIC (Applicata direttamente alla cella)
%         if ~isnan(subsectionData.Telaps(s))
%             % LEVEL 1: CLASSIC BPT (Data conosciuta, es. Vettore = 10, Fucino = 111)
%             q1 = funBPT(subsectionData.Telaps(s));
%             q2 = funBPT(subsectionData.Telaps(s) + t);
%             if q1 < 0.9999
%                 subBPTProb(s) = (q2 - q1) / (1 - q1);
%             else
%                 subBPTProb(s) = 1.0;
%             end
% 
%         elseif ~isnan(subsectionData.Thist(s))
%             % LEVEL 2: HISTORICAL OPEN INTERVAL
%             q1 = integral(funBPT, subsectionData.Thist(s), subsectionData.Thist(s) + t);
%             funDenom = @(T) (1 - cdf('inversegaussian', T, mu_s, lambda_s));
%             q2 = integral(funDenom, subsectionData.Thist(s), inf);      
%             subBPTProb(s) = (t - q1) / q2;
% 
%         else
%             % LEVEL 3: UNKNOWN DATE
%             q_unk = integral(funBPT, 0, t);                      
%             subBPTProb(s) = (t - q_unk) / mu_s;
%         end
% 
%         % Safety Checks
%         if subBPTProb(s) > 1, subBPTProb(s) = 1; end
%         if subBPTProb(s) < 0 || isnan(subBPTProb(s)), subBPTProb(s) = 0; end
%     end
%     warning('on', id_warn);
% 
%     % --- Step 3: Aggregate Parent Fault Probabilities ---
%     faultParticipationRate = zeros(numFaults, 1);
%     faultPoissonProb       = zeros(numFaults, 1);
%     faultBPTProb           = zeros(numFaults, 1);
% 
%     % Build the global lookup map
%     globalLookup = struct('nf', cell(numSubsections, 1));
%     current_global = 1;
%     for nf = 1:numFaults
%         for s_idx = 1:faultResults(nf).m
%             for d_idx = 1:faultResults(nf).n
%                 globalLookup(current_global).nf = nf;
%                 current_global = current_global + 1;
%             end
%         end
%     end
% 
%     for nf = 1:numFaults
%         subIdx = find([globalLookup.nf] == nf);
%         if isempty(subIdx), continue; end
% 
%         % Calculate total fault rate for events >= mw_lim
%         for r = 1:length(bestRates.mean)
%             if bestRates.mean(r) > 0 && invData.Mw(r) >= mw_lim
%                 if any(logicalGsr(r, subIdx))
%                     faultParticipationRate(nf) = faultParticipationRate(nf) + bestRates.mean(r);
%                 end
%             end
%         end
% 
%         % Parent Fault Poisson Probability
%         faultPoissonProb(nf) = 1 - exp(-t * faultParticipationRate(nf));
% 
%         % Parent Fault BPT Aggregation via Subsection Union Rule
%         % Uniamo le probabilità condizionate delle sottosezioni che compongono la faglia
%         faultBPTProb(nf) = 1 - prod(1 - subBPTProb(subIdx));
%     end
% 
%     % --- Step 4: Compile Outputs Tables ---
%     outputs = struct();
%     outputs.Subsections = table((1:numSubsections)', subEventRate(:), subTmean(:), ...
%         subPoissonProb(:), subBPTProb(:), ...
%         'VariableNames', {'Subsection_ID', 'ParticipationRate', 'Tmean', 'PoissonProb', 'BPT_TimeDependentProb'});
% 
%     outputs.Faults = table({faultResults.Name}', faultParticipationRate, faultPoissonProb, faultBPTProb, ...
%         'VariableNames', {'FaultName', 'ParticipationRate', 'PoissonProb', 'BPT_TimeDependentProb'});
% end