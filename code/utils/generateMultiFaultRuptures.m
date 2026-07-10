function All_Final_Ruptures = generateMultiFaultRuptures(All_Elementary_Ruptures, Global_Subsections_Data, A_S2S, params)
% Recursively generates complex multi-fault ruptures.
% Uses a Core-Tapering expansion approach restricted by structural connectivity
% and down-dip alignment constraints. 

faultIDs = sort(unique(Global_Subsections_Data.Fault_ID));
numFaults = max(faultIDs);

% --- 1. PRE-COMPUTE F2F CONNECTIVITY ---
F2F_Conn = cell(numFaults, numFaults);
for i = 1:numel(faultIDs)
    for j = i+1:numel(faultIDs)
        fA = faultIDs(i); fB = faultIDs(j);
        idxA = find(Global_Subsections_Data.Fault_ID == fA);
        idxB = find(Global_Subsections_Data.Fault_ID == fB);
        [subA, subB] = find(A_S2S(idxA, idxB));
        if ~isempty(subA)
            F2F_Conn{fA, fB} = [idxA(subA), idxB(subB)];
            F2F_Conn{fB, fA} = [idxB(subB), idxA(subA)];
        end
    end
end

% --- 2. PRE-COMPUTE FULL-DIP RUPTURES ---
FullDipRuptures = cell(numFaults, 1);
for f = 1:numel(faultIDs)
    fid = faultIDs(f);
    R_f = All_Elementary_Ruptures(All_Elementary_Ruptures.Fault_ID == fid, :);
    [~, maxIdx] = max(R_f.Area); 
    FullDipRuptures{fid} = R_f(maxIdx, :);
end

% --- 3. RECURSIVE SEARCH WITH INTEGRATED FILTERS ---
max_est = 2000000;
multi_results = cell(max_est, 5);
rupt_count = 0;

fprintf('🚀 Starting Optimized Core-Tapering Expansion...\n');

for f = 1:numel(faultIDs)
    if mod(f, 5) == 0 || f == 1
        fprintf('   ...processing fault %d of %d (Current count: %d)\n', f, numel(faultIDs), rupt_count);
    end
    start_fid = faultIDs(f);
    R_start = All_Elementary_Ruptures(All_Elementary_Ruptures.Fault_ID == start_fid, :);
    
    for r = 1:height(R_start)
        [multi_results, rupt_count] = expandCoreRecursive(...
            start_fid, start_fid, R_start.Global_IDs{r}, R_start.Area(r), R_start.Rupture_ID(r), ...
            F2F_Conn, All_Elementary_Ruptures, FullDipRuptures, params, Global_Subsections_Data, ...
            multi_results, rupt_count);
    end
end

% --- 4. FINALIZE & UNIQUE ---
if rupt_count > 0
    Raw_Table = cell2table(multi_results(1:rupt_count, :), ...
        'VariableNames', {'Rupture_ID', 'Elementary_IDs', 'Area', 'Mw', 'Global_IDs'});
    
   % Ensure Elementary_IDs is a cell array (even for cases with only 2 faults connected)
if ~iscell(Raw_Table.Elementary_IDs)
    e_ids = num2cell(Raw_Table.Elementary_IDs, 2);
else
    e_ids = Raw_Table.Elementary_IDs;
end

% Sort IDs internally to prepare for unique pairing checks
sorted_EIDs = cellfun(@(x) sort(x), e_ids, 'UniformOutput', false);
    unique_keys = string(cellfun(@(x) strjoin(string(x), ','), sorted_EIDs, 'UniformOutput', false));
    [~, uIdx] = unique(unique_keys, 'stable');
    
    All_Final_Ruptures = Raw_Table(uIdx, :);
    All_Final_Ruptures.Rupture_ID = (1:height(All_Final_Ruptures))';
    
    fprintf('✅ Done. Generated %d unique ruptures (Mw filtered, Dip-aligned).\n', height(All_Final_Ruptures));
else
    All_Final_Ruptures = table();
    fprintf('⚠️ No ruptures generated.\n');
end
end

%% ================== OPTIMIZED RECURSIVE FUNCTION ==================
function [list, count] = expandCoreRecursive(last_fid, path_fids, gIDs, area, eIDs, ...
    F2F, All_Elem, FullDip, params, GSD, list, count)

% PRUNING Mw: If the candidate exceeds maxMw, terminate further branch exploration.
current_mw = seismicScalingEngine(area, 'Area2Mw', params.faultMechanism);
if current_mw > params.maxMw, return; end

neighbors = find(~cellfun(@isempty, F2F(last_fid, :)));
neighbors = setdiff(neighbors, path_fids);
 
% Pre-extract current Dip levels for downward alignment filtering
current_fault_gIDs = gIDs(ismember(gIDs, find(GSD.Fault_ID == last_fid)));
dips_current = unique(GSD.D_idx(ismember(GSD.Global_ID, current_fault_gIDs)));

for n = 1:numel(neighbors)
    neighbor_fid = neighbors(n);
    bridge_data = F2F{last_fid, neighbor_fid};
    
    if isempty(intersect(gIDs, bridge_data(:,1))), continue; end
    
    % OPTION A: TERMINAL (Tapered)
    R_neighbor = All_Elem(All_Elem.Fault_ID == neighbor_fid, :);
    for rn = 1:height(R_neighbor)
        if any(ismember(R_neighbor.Global_IDs{rn}, bridge_data(:,2)))
            
            % Check Dip Alignment
            dips_neighbor = unique(GSD.D_idx(ismember(GSD.Global_ID, R_neighbor.Global_IDs{rn})));
            
            if ~isempty(intersect(dips_current, dips_neighbor))
                new_area = area + R_neighbor.Area(rn);
                new_mw = seismicScalingEngine(new_area, 'Area2Mw', params.faultMechanism);
                
                if new_mw <= params.maxMw
                    count = count + 1;
                    if count > size(list,1), list = [list; cell(50000, 5)]; end
                    
                    new_eIDs = [eIDs, R_neighbor.Rupture_ID(rn)];
                    new_gIDs = unique([gIDs; R_neighbor.Global_IDs{rn}]);
                    list(count, :) = {count, new_eIDs, new_area, new_mw, new_gIDs};
                end
            end
        end
    end
    
    % OPTION B: INTERNAL (Full-Dip)
    R_full = FullDip{neighbor_fid};
    if any(ismember(R_full.Global_IDs{1}, bridge_data(:,2)))
        new_area_full = area + R_full.Area(1);
        new_eIDs_full = [eIDs, R_full.Rupture_ID(1)];
        new_gIDs_full = unique([gIDs; R_full.Global_IDs{1}]);
        
        [list, count] = expandCoreRecursive(...
            neighbor_fid, [path_fids, neighbor_fid], new_gIDs_full, new_area_full, ...
            new_eIDs_full, F2F, All_Elem, FullDip, params, GSD, list, count);
    end
end
end