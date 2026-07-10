function A_S2S = computeS2SMatrix(faultResults, Global_Subsections_Data, params)
% COMPUTES2SMATRIX - Computes the Subsection-to-Subsection connectivity matrix.
% This matrix defines "bridges" where a rupture can propagate between different faults.
%
% Optimization: Uses SystemID grouping and spatial filtering to avoid O(n^2) complexity.

    fprintf('🔗 Building Subsection-to-Subsection (S2S) connectivity matrix...\n');

    % --- Parameters ---
    MAX_GAP_SHALLOW = params.maxGapShallow;
    MAX_GAP_DEEP    = params.maxGapDeep;
    Z_LIMIT         = params.zLimit;
    
    nTotalSubs = size(Global_Subsections_Data, 1);
    A_S2S = sparse(nTotalSubs, nTotalSubs);
    
    % Group faults by their SystemID
    systemIDs = [faultResults.SystemID];
    uniqueSystems = unique(systemIDs);
    
    % Indexing helper: Map Fault IDs to Global Indices
    faultIDCol = Global_Subsections_Data.Fault_ID;

    for s = 1:length(uniqueSystems)
        sysID = uniqueSystems(s);
        faultsInSys = find(systemIDs == sysID);
        
        if length(faultsInSys) < 2, continue; end % Skip isolated faults
        
        % Iterate over pairs within the same system
        for i = 1:length(faultsInSys)
            for j = i+1:length(faultsInSys)
                fidA = faultsInSys(i);
                fidB = faultsInSys(j);
                
                % Get subsection data for both faults
                idxA = find(faultIDCol == fidA);
                idxB = find(faultIDCol == fidB);
                subA = Global_Subsections_Data(idxA, :);
                subB = Global_Subsections_Data(idxB, :);
                
                % --- STEP 1: FIND NEAREST STRIKE COLUMNS ---
                % We identify which column of Fault A is closest to which column of Fault B
                [bestSA, bestSB] = findClosestColumns(subA, subB);
                
                if isnan(bestSA), continue; end
                
                % --- STEP 2: CONNECT SUBSECTIONS WITHIN BEST COLUMNS ---
                % Filter subsets to only these columns
                colA_logic = (subA.S_idx == bestSA);
                colB_logic = (subB.S_idx == bestSB);
                
                colA = subA(colA_logic, :);
                colB = subB(colB_logic, :);
                
                % Global indices for the sparse matrix
                globalIdxA = idxA(colA_logic);
                globalIdxB = idxB(colB_logic);
                
                n_conn = 0;
                for rowA = 1:height(colA)
                    for rowB = 1:height(colB)
                        
                        % Rule: Only connect same depth level (D_idx)
                        if colA.D_idx(rowA) ~= colB.D_idx(rowB), continue; end
                        
                        % Distance calculation between top edges
                        d = min([
                            sqrt((colA.X1_km(rowA)-colB.X1_km(rowB))^2 + (colA.Y1_km(rowA)-colB.Y1_km(rowB))^2);
                            sqrt((colA.X1_km(rowA)-colB.X2_km(rowB))^2 + (colA.Y1_km(rowA)-colB.Y2_km(rowB))^2);
                            sqrt((colA.X2_km(rowA)-colB.X1_km(rowB))^2 + (colA.Y2_km(rowA)-colB.Y1_km(rowB))^2);
                            sqrt((colA.X2_km(rowA)-colB.X2_km(rowB))^2 + (colA.Y2_km(rowA)-colB.Y2_km(rowB))^2)
                        ]);
                        
                        % Variable Threshold based on depth
                        zRef = colA.Z_km(rowA);
                        maxGap = ifthen(zRef <= Z_LIMIT, MAX_GAP_SHALLOW, MAX_GAP_DEEP);
                        
                        if d <= maxGap
                            A_S2S(globalIdxA(rowA), globalIdxB(rowB)) = 1;
                            A_S2S(globalIdxB(rowB), globalIdxA(rowA)) = 1;
                            n_conn = n_conn + 1;
                        end
                    end
                end
                if n_conn > 0
                    fprintf('   🔗 Connected: %s <-> %s (%d bridges)\n', ...
                        faultResults(fidA).Name, faultResults(fidB).Name, n_conn);
                end
            end
        end
    end
    fprintf('✅ S2S Matrix built. Total connections: %d\n', nnz(A_S2S)/2);
end

% --- Helper Functions ---

function [bestSA, bestSB] = findClosestColumns(subA, subB)
    % Find the pair of strike indices (S_idx) with the minimum 3D distance
    uSA = unique(subA.S_idx);
    uSB = unique(subB.S_idx);
    minD = inf;
    bestSA = NaN; bestSB = NaN;
    
    for sA = uSA'
        dataA = subA(subA.S_idx == sA, :);
        for sB = uSB'
            dataB = subB(subB.S_idx == sB, :);
            
            % Check first depth level of each column for speed
            d = sqrt((dataA.X_km(1)-dataB.X_km(1))^2 + (dataA.Y_km(1)-dataB.Y_km(1))^2);
            if d < minD
                minD = d;
                bestSA = sA;
                bestSB = sB;
            end
        end
    end
end

function val = ifthen(condition, trueVal, falseVal)
    if condition, val = trueVal; else, val = falseVal; end
end