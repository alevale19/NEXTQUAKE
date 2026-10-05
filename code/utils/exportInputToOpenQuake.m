function exportInputToOpenQuake(bestRates, outputs, invData, All_Rupture_Candidates, faultResults, params)
% EXPORTINPUTTOOPENQUAKE - Generates OpenQuake source files (Poisson always, BPT if available)
% by reading exact geographic mesh vertices from faultResults and output rates.

    fprintf('🌐 Exporting active ruptures to OpenQuake SourceModel format...\n');
    
    % --- Step A: Directory Setup ---
    outputDir = './Output/OpenQuake_Ruptures';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    
    % Rate threshold for active ruptures
    rateThreshold = 1e-8; 
    
    % Extract Poisson rates
    poissonRates = bestRates.mean; 
    
    % Check if BPT rates are available (case-insensitive check)
    hasBPT = false;
    bptRates = [];
    
    if isstruct(outputs)
        fields = fieldnames(outputs);
        matchIdx = find(strcmpi(fields, 'bestRates_BPT'), 1);
        if ~isempty(matchIdx)
            hasBPT = true;
            bptRates = outputs.(fields{matchIdx});
            fprintf('  -> BPT rates found! Generating both Poisson and BPT models.\n');
        end
    end
    
    if ~hasBPT
        fprintf('  -> BPT rates NOT found in outputs. Skipping BPT export (Poisson only).\n');
    end
    
    % Identify active ruptures
    if hasBPT
        activeRupIdx = find(poissonRates > rateThreshold | bptRates > rateThreshold);
    else
        activeRupIdx = find(poissonRates > rateThreshold);
    end
    
    numActive = numel(activeRupIdx);
    fprintf('  -> Found %d significant active ruptures (rate > %.1e).\n', ...
            numActive, rateThreshold);
        
    % --- Step B: Pre-index Mapping (Global ID -> faultIndex, s, d) ---
    globalLookup = struct('nf', {}, 's', {}, 'd', {});
    current_global = 1;
    for nf = 1:length(faultResults)
        for s = 1:faultResults(nf).m
            for d = 1:faultResults(nf).n
                globalLookup(current_global).nf = nf;
                globalLookup(current_global).s  = s;
                globalLookup(current_global).d  = d;
                current_global = current_global + 1;
            end
        end
    end
    
    % --- Step C: XML Document Initialization ---
    [docNode_P, srcGroupNode_P] = createBaseOQDocument('Poisson Model');
    
    if hasBPT
        [docNode_BPT, srcGroupNode_BPT] = createBaseOQDocument('BPT Time-Dependent Model');
    end
    
    % --- Step D: Loop over active ruptures ---
    for i = 1:numActive
        rupID = activeRupIdx(i);
        rateP = poissonRates(rupID);
        mw    = invData.Mw(rupID);
        
        if hasBPT
            rateBPT = bptRates(rupID);
        end
        
        % Extract indices of involved subsections
        if iscell(All_Rupture_Candidates)
            subSecIndices = All_Rupture_Candidates{rupID}; 
        elseif istable(All_Rupture_Candidates)
            if ismember('Global_IDs', All_Rupture_Candidates.Properties.VariableNames)
                subSecIndices = All_Rupture_Candidates.Global_IDs{rupID};
            else
                subSecIndices = All_Rupture_Candidates.global_ids{rupID};
            end
        else
            subSecIndices = All_Rupture_Candidates(rupID, :);
        end
        
        if iscell(subSecIndices), subSecIndices = [subSecIndices{:}]; end
        subSecIndices = double(subSecIndices(:)');
        
        posListStr = ''; 
        totalArea = 0;
        weightedRakeSum = 0;
        
        % Extract geometry for the rupture
        for s_idx_map = 1:numel(subSecIndices)
            gID = subSecIndices(s_idx_map);
            
            % Retrieve fault index and matrix grid indices (s, d)
            nf = globalLookup(gID).nf;
            s  = globalLookup(gID).s;
            d  = globalLookup(gID).d;
            
            f = faultResults(nf);
            
            % Extract meshes
            lat_m   = f.lat_mesh;
            lon_m   = f.lon_mesh;
            depth_m = f.depth_mesh;
            
            % Exact 4 Corner Mesh Extractor
            Lon_v = [lon_m(s, d), lon_m(s+1, d), lon_m(s+1, d+1), lon_m(s, d+1)];
            Lat_v = [lat_m(s, d), lat_m(s+1, d), lat_m(s+1, d+1), lat_m(s, d+1)];
            Z_v   = [depth_m(s, d), depth_m(s+1, d), depth_m(s+1, d+1), depth_m(s, d+1)];
            
            % Append vertices to OpenQuake posList format ("Lon Lat Z ")
            for v = 1:4
                posListStr = [posListStr, sprintf('%.5f %.5f %.2f ', Lon_v(v), Lat_v(v), Z_v(v))]; %#ok<AGROW>
            end
            
            % Calculate area-weighted mean Rake
            cell_len = hypot(f.coordsUTM(min(s+1,end),1) - f.coordsUTM(s,1), ...
                             f.coordsUTM(min(s+1,end),2) - f.coordsUTM(s,2)) / 1000; 
            cell_w = f.Width / f.n;
            area_s = cell_len * cell_w;
            
            weightedRakeSum = weightedRakeSum + (f.Rake * area_s);
            totalArea = totalArea + area_s;
        end
        
        meanRake = weightedRakeSum / totalArea;
        posListStr = strtrim(posListStr); 
        
        % Append Rupture Node to Poisson File
        if rateP > rateThreshold
            appendRuptureNode(docNode_P, srcGroupNode_P, rupID, rateP, mw, meanRake, posListStr);
        end
        
        % Append Rupture Node to BPT File (only if available and active)
        if hasBPT && rateBPT > rateThreshold
            appendRuptureNode(docNode_BPT, srcGroupNode_BPT, rupID, rateBPT, mw, meanRake, posListStr);
        end
    end
    
    % --- Step E: Save Final Files ---
    filePoisson = fullfile(outputDir, 'OQ_Source_Model_Poisson.xml');
    xmlwrite(filePoisson, docNode_P);
    fprintf('✅ Poisson OpenQuake XML model successfully saved at: %s\n', filePoisson);
    
    if hasBPT
        fileBPT = fullfile(outputDir, 'OQ_Source_Model_BPT.xml');
        xmlwrite(fileBPT, docNode_BPT);
        fprintf('✅ BPT OpenQuake XML model successfully saved at:     %s\n', fileBPT);
    end
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [docNode, srcGroupNode] = createBaseOQDocument(modelName)
    docNode = com.mathworks.xml.XMLUtils.createDocument('nrml');
    nrmlNode = docNode.getDocumentElement;
    nrmlNode.setAttribute('xmlns', 'http://openquake.org/xmlns/nrml/0.5');
    nrmlNode.setAttribute('xmlns:gml', 'http://www.opengis.net/gml');
    
    srcModelNode = docNode.createElement('sourceModel');
    srcModelNode.setAttribute('name', modelName);
    nrmlNode.appendChild(srcModelNode);
    
    srcGroupNode = docNode.createElement('sourceGroup');
    srcGroupNode.setAttribute('name', 'rupture group');
    srcGroupNode.setAttribute('tectonicRegion', 'Active Shallow Crust');
    srcModelNode.appendChild(srcGroupNode);
end

function appendRuptureNode(docNode, srcGroupNode, rupID, rate, mw, meanRake, posListStr)
    faultSrcNode = docNode.createElement('characteristicFaultSource');
    faultSrcNode.setAttribute('id', sprintf('%d', rupID));
    faultSrcNode.setAttribute('name', sprintf('rupture %d', rupID));
    srcGroupNode.appendChild(faultSrcNode);
    
    mfdNode = docNode.createElement('arbitraryMFD');
    faultSrcNode.appendChild(mfdNode);
    
    occurRatesNode = docNode.createElement('occurRates');
    occurRatesNode.appendChild(docNode.createTextNode(sprintf('%.2e', rate)));
    mfdNode.appendChild(occurRatesNode);
    
    magnitudesNode = docNode.createElement('magnitudes');
    magnitudesNode.appendChild(docNode.createTextNode(sprintf('%.3f', mw)));
    mfdNode.appendChild(magnitudesNode);
    
    rakeNode = docNode.createElement('rake');
    rakeNode.appendChild(docNode.createTextNode(sprintf('%.1f', meanRake)));
    faultSrcNode.appendChild(rakeNode);
    
    surfaceWrapperNode = docNode.createElement('surface');
    faultSrcNode.appendChild(surfaceWrapperNode);
    
    gridSurfNode = docNode.createElement('griddedSurface');
    surfaceWrapperNode.appendChild(gridSurfNode);
    
    posListNode = docNode.createElement('gml:posList');
    posListNode.appendChild(docNode.createTextNode(posListStr));
    gridSurfNode.appendChild(posListNode);
end
