function output = seismicScalingEngine(inputVal, mode, faultMechanism)
% Inputs:
%    inputVal       - Numeric input value:
%                     -> Magnitude (Mw) if mode is 'Mw2Length'
%                     -> Rupture area (km²) if mode is 'Area2Mw'
%    mode           - Control string:
%                     -> 'Mw2Length' : Computes target subsection length (km)
%                     -> 'Area2Mw'   : Computes consistent Mw from Area (km²)
%    faultMechanism - Fault kinematic style: 
%                     -> 'All' (Crustal Ensemble), 'N' (Normal), 'SS' (Strike-Slip), 
%                        'R' (Reverse), 'SUB' (Subduction Interface)

    % Handle optional fault mechanism parameter (defaults to 'All')
    if nargin < 3, faultMechanism = 'All'; end

    switch mode
        
        case 'Mw2Length'
            % =============================================================
            % DIRECT MODE: Magnitude (Mw) -> Subsection Length (km)
            % =============================================================
            % Computes theoretical Area and extracts its square root (equivalent length)
            
            % 1. Wells & Coppersmith (1994)
            %L_RA_WC94All = sqrt(10^(-3.49 + 0.91 * inputVal)); % Global crustal model
            L_RA_WC94N   = sqrt(10^(-2.87 + 0.82 * inputVal)); % Specific for Normal (extensional) faults
            L_RA_WC94SS  = sqrt(10^(-3.42 + 0.90 * inputVal)); % Specific for Strike-Slip faults
            L_RA_WC94R   = sqrt(10^(-3.99 + 0.98 * inputVal)); % Specific for Reverse faults

            % 2. Thingbaijam et al. (2017)
            L_RA_T17sub  = sqrt(10^(-3.292 + 0.949 * inputVal)); % Specific for Subduction zones
            L_RA_T174N   = sqrt(10^(-2.551 + 0.808 * inputVal)); % Specific for Normal faults
            L_RA_T17SS   = sqrt(10^(-3.486 + 0.942 * inputVal)); % Specific for Strike-Slip faults
            L_RA_T17R    = sqrt(10^(-4.362 + 1.049 * inputVal)); % Specific for Reverse (shallow crustal) faults

            % Filter and ensemble models based on fault kinematics
            switch upper(faultMechanism)
                case 'N'     % Only Normal faulting relations
                    selected_models = [L_RA_WC94N, L_RA_T174N];
                case 'SS'    % Only Strike-Slip relations
                    selected_models = [L_RA_WC94SS, L_RA_T17SS];
                case 'R'     % Only Reverse / Thrust relations
                    selected_models = [L_RA_WC94R, L_RA_T17R];
                case 'SUB'   % Subduction Zones: strictly Thingbaijam 2017 orthogonal relation
                    selected_models = L_RA_T17sub;
                otherwise    % Crustal ensemble (Subduction is excluded here)
                    selected_models = [L_RA_WC94N, L_RA_WC94SS, L_RA_WC94R, ...
                                       L_RA_T174N, L_RA_T17SS, L_RA_T17R];
            end
            
            % Compute the arithmetic mean of the ensemble (output is in kilometers)
            output = mean(selected_models);

            fprintf('=== Scaling Engine Report: Mw -> Length ===\n');
            fprintf('Input Forecast Mw:           %.2f\n', inputVal);
            fprintf('Selected Mechanism:          %s\n', faultMechanism);
            fprintf('Computed Subsection Length:  %.2f km\n', output);
            fprintf('===========================================\n\n');

        case 'Area2Mw'
            % =============================================================
            % INVERSE MODE: Rupture Area (km²) -> Magnitude (Mw)
            % =============================================================
            % Applies native inverse regressions originally fit by the authors on catalogs
            logA = log10(inputVal);
            
            % 1. Wells & Coppersmith (1994)
            %Mw_WC94All = 4.07 + 0.98 * logA; % Global crustal inverse model
            Mw_WC94N   = 3.93 + 1.02 * logA; % Specific for Normal faults
            Mw_WC94SS  = 3.98 + 1.02 * logA; % Specific for Strike-Slip faults
            Mw_WC94R   = 4.33 + 0.90 * logA; % Specific for Reverse faults

            % 2. Thingbaijam et al. (2017)
            Mw_T174N   = (logA - (-2.551)) / 0.808; % Specific for Normal faults
            Mw_T17SS   = (logA - (-3.486)) / 0.942; % Specific for Strike-Slip faults
            Mw_T17R    = (logA - (-4.362)) / 1.049; % Specific for Reverse (shallow crustal) faults
            Mw_T17sub  = (logA - (-3.292)) / 0.949; % Specific for Subduction zones

            % Filter and ensemble corresponding inverse equations
            switch upper(faultMechanism)
                case 'N'     % Synchronized with the Apenninic sampling strategy
                    selected_models = [Mw_WC94N, Mw_T174N];
                case 'SS'
                    selected_models = [Mw_WC94SS, Mw_T17SS];
                case 'R'
                    selected_models = [Mw_WC94R, Mw_T17R];
                case 'SUB'   % Subduction Zones: strictly Thingbaijam 2017 orthogonal relation
                    selected_models = Mw_T17sub;
                otherwise    % Crustal inverse ensemble (Subduction is excluded here)
                    selected_models = [Mw_WC94N, Mw_WC94SS, Mw_WC94R, ...
                                       Mw_T174N, Mw_T17SS, Mw_T17R];
            end
            
            % Compute the ensemble average Magnitude
            output = mean(selected_models);

        otherwise
            % Safety block for potential typos in the main codebase calling routine
            error('Unrecognized mode! Choose strictly between ''Mw2Length'' or ''Area2Mw''.');
    end
end