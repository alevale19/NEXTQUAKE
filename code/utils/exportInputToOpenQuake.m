function exportInputToOpenQuake(bestRates, invData, All_Rupture_Candidates, faultResults, params)
% EXPORTALLRUPTURESTOOPENQUAKE - Generates a single unified OpenQuake source file 
% by reading the exact 4 geographic mesh vertices (lat_mesh, lon_mesh) for each 
% subsection directly from the faultResults structure array.

    fprintf('🌐 Exporting active ruptures to OpenQuake SourceModel format...\n');
    
    % --- Step A: Directory and Filename Setup ---
    outputDir = './Output/OpenQuake_Ruptures';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    outputFile = fullfile(outputDir, 'OQ_Source_Model.xml');
    
    % Threshold to filter out virtually zero-activity ruptures
    rateThreshold = 1e-8; 
    activeRupIdx = find(bestRates.mean > rateThreshold);
    numActive = numel(activeRupIdx);
    
    fprintf('  -> Found %d significant active ruptures (rate > %.1e).\n', ...
            numActive, rateThreshold);
        
    % --- Step B: XML Document Initialization ---
    docNode = com.mathworks.xml.XMLUtils.createDocument('nrml');
    nrmlNode = docNode.getDocumentElement;
    nrmlNode.setAttribute('xmlns', 'http://openquake.org/xmlns/nrml/0.5');
    nrmlNode.setAttribute('xmlns:gml', 'http://www.opengis.net/gml');
    
    srcModelNode = docNode.createElement('sourceModel');
    srcModelNode.setAttribute('name', 'Converted Rupture Model');
    nrmlNode.appendChild(srcModelNode);
    
    srcGroupNode = docNode.createElement('sourceGroup');
    srcGroupNode.setAttribute('name', 'rupture group');
    srcGroupNode.setAttribute('tectonicRegion', 'Active Shallow Crust');
    srcModelNode.appendChild(srcGroupNode);
    
    % --- Step C: Loop over active ruptures ---
    for i = 1:numActive
        rupID = activeRupIdx(i);
        rate = bestRates.mean(rupID);
        mw = invData.Mw(rupID);
        
        % Extract indices of the involved subsections
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
        
        % --- Pre-index Mapping ---
        % Since subSecIndices references a global incremental ID (1 to total_subs),
        % we build a quick structural lookup map matching global_id -> [fault_index, s, d]
        % exactly matching the generation loop order in computeSubsectionCentroids.
        globalLookup = [];
        current_global = 1;
        for nf = 1:length(faultResults)
            for s = 1:faultResults(nf).m
                for d = 1:faultResults(nf).n
                    globalLookup(current_global).nf = nf;
                    globalLookup(current_global).s = s;
                    globalLookup(current_global).d = d;
                    current_global = current_global + 1;
                end
            end
        end
        
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
            
            % =============================================================
            % EXACT 4 CORNER MESH EXTRACTOR (Lon, Lat, Depth)
            % Matches the exact topology of the subsection discretization grid
            % =============================================================
            % Corner 1: Top Left
            Lon_v(1) = lon_m(s, d);     Lat_v(1) = lat_m(s, d);     Z_v(1) = depth_m(s, d);
            % Corner 2: Top Right
            Lon_v(2) = lon_m(s+1, d);   Lat_v(2) = lat_m(s+1, d);   Z_v(2) = depth_m(s+1, d);
            % Corner 3: Bottom Right
            Lon_v(3) = lon_m(s+1, d+1); Lat_v(3) = lat_m(s+1, d+1); Z_v(3) = depth_m(s+1, d+1);
            % Corner 4: Bottom Left
            Lon_v(4) = lon_m(s, d+1);   Lat_v(4) = lat_m(s, d+1);   Z_v(4) = depth_m(s, d+1);
            
            % Append vertices to OpenQuake posList format ("Lon Lat Z ")
            for v = 1:4
                posListStr = [posListStr, sprintf('%.5f %.5f %.2f ', Lon_v(v), Lat_v(v), Z_v(v))];
            end
            
            % Calculate area-weighted mean Rake
            % Area = cell_length * cell_width
            cell_len = hypot(f.coordsUTM(min(s+1,end),1) - f.coordsUTM(s,1), ...
                             f.coordsUTM(min(s+1,end),2) - f.coordsUTM(s,2)) / 1000; 
            cell_w = f.Width / f.n;
            area_s = cell_len * cell_w;
            
            rake_s = f.Rake;
            weightedRakeSum = weightedRakeSum + (rake_s * area_s);
            totalArea = totalArea + area_s;
        end
        
        meanRake = weightedRakeSum / totalArea;
        posListStr = strtrim(posListStr); 
        
        % 4. Build XML structured nodes
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
    
    % --- Step D: Save Final File ---
    xmlwrite(outputFile, docNode);
    fprintf('✅ Unified OpenQuake XML source model successfully generated at: %s\n', outputFile);
end
