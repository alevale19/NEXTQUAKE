% =========================================================================
% SEISMIC SCALING RELATIONS & SUBSECTION SIZING ANALYSIS
% =========================================================================

clear; clc; close all;

%% 1. PARAMETERS & INPUT MAGNITUDE RANGE
forecastMw = 4.5:0.1:9.2; % Magnitude range for evaluation
inputVal   = forecastMw;  % Alias for consistency

%% 2. EMPIRICAL SCALING RELATIONS (RUPTURE AREA - RA [km^2])
% Wells & Coppersmith (1994) - Crustal faults
L_RA_WC94N  = 10.^(-2.87 + 0.82 * inputVal); % Normal faults (Mw 5.2 - 7.3)
L_RA_WC94SS = 10.^(-3.42 + 0.90 * inputVal); % Strike-Slip faults (Mw 4.8 - 7.9)
L_RA_WC94R  = 10.^(-3.99 + 0.98 * inputVal); % Reverse faults (Mw 4.8 - 7.6)

% Thingbaijam et al. (2017) - Crustal & Subduction faults
L_RA_T17N   = 10.^(-2.551 + 0.808 * inputVal); % Normal faults (Mw 5.8 - 8.4)
L_RA_T17SS  = 10.^(-3.486 + 0.942 * inputVal); % Strike-Slip faults (Mw 5.4 - 8.7)
L_RA_T17R   = 10.^(-4.362 + 1.049 * inputVal); % Reverse faults (Mw 5.6 - 7.7)
L_RA_T17sub = 10.^(-3.292 + 0.949 * inputVal); % Subduction zones (Mw 6.7 - 9.2)

%% 3. ENSEMBLE MEAN CHARACTERISTIC LENGTHS (sqrt(RA) [km])
av_N   = (sqrt(L_RA_WC94N)  + sqrt(L_RA_T17N))  / 2;
av_R   = (sqrt(L_RA_WC94R)  + sqrt(L_RA_T17R))  / 2;
av_SS  = (sqrt(L_RA_WC94SS) + sqrt(L_RA_T17SS)) / 2;
av_ALL = (sqrt(L_RA_WC94N)  + sqrt(L_RA_T17N) + ...
          sqrt(L_RA_WC94R)  + sqrt(L_RA_T17R) + ...
          sqrt(L_RA_WC94SS) + sqrt(L_RA_T17SS)) / 6;

%% 4. SEISMIC MOMENT & AVERAGE DISPLACEMENT CALCULATION
Mo       = 10.^(1.5 * inputVal + 9.05); % Seismic Moment [N*m]
areas_m2 = L_RA_T17SS * 1e6;            % Convert km^2 to m^2 (using T17 SS Area)
mu       = 3e10;                        % Shear modulus [Pa]
Dr       = Mo ./ (mu * areas_m2);       % Average Displacement [m]

%% 5. MULTI-PANEL VISUALIZATION
figure('Name', 'Scaling Relations Calibration & Subsection Sizing', 'Color', 'w', 'Position', [100 100 1200 450]);

% --- Subplot 1: Rupture Area vs. Magnitude (Validity Ranges) ---
subplot(1, 3, 1);
% Wells & Coppersmith (1994) - Solid lines
semilogx(L_RA_WC94N(8:29),   inputVal(8:29),   'k',  'LineWidth', 2, 'DisplayName', 'WC94 Normal'); hold on;
semilogx(L_RA_WC94SS(4:35),  inputVal(4:35),  'r',  'LineWidth', 2, 'DisplayName', 'WC94 SS');
semilogx(L_RA_WC94R(4:32),   inputVal(4:32),   'b',  'LineWidth', 2, 'DisplayName', 'WC94 Reverse');

% Thingbaijam et al. (2017) - Dashed lines
semilogx(L_RA_T17N(14:36),   inputVal(14:36),  'k--', 'LineWidth', 2, 'DisplayName', 'T17 Normal');
semilogx(L_RA_T17SS(10:36),  inputVal(10:36),  'r--', 'LineWidth', 2, 'DisplayName', 'T17 SS');
semilogx(L_RA_T17R(12:32),   inputVal(12:32),  'b--', 'LineWidth', 2, 'DisplayName', 'T17 Reverse');
semilogx(L_RA_T17sub(23:48), inputVal(23:48), 'm--', 'LineWidth', 2, 'DisplayName', 'T17 Subduction');

axis square; grid on; box on;
ylabel('M_w', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Area (km^2)', 'FontSize', 11, 'FontWeight', 'bold');
ylim([4.5 8.5]);
title('a) Empirical Rupture Area', 'FontSize', 12);

% --- Subplot 2: Average Displacement vs. Magnitude ---
subplot(1, 3, 2);
plot(Dr, inputVal, 'k', 'LineWidth', 2);
axis square; grid on; box on;
ylabel('M_w', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('D (m)', 'FontSize', 11, 'FontWeight', 'bold');
xticks(0:2:16);
ylim([4.5 8.5]);
title('b) Average Slip (D)', 'FontSize', 12);

% --- Subplot 3: Ensemble Mean Characteristic Length ---
subplot(1, 3, 3);
plot(av_N,                 forecastMw, 'r',  'LineWidth', 2, 'DisplayName', 'Normal'); hold on;
plot(av_SS,                forecastMw, 'b',  'LineWidth', 2, 'DisplayName', 'Strike-Slip');
plot(av_R,                 forecastMw, 'g',  'LineWidth', 2, 'DisplayName', 'Reverse');
plot(av_ALL,               forecastMw, 'r:', 'LineWidth', 2, 'DisplayName', 'All Crustal');
plot(sqrt(L_RA_T17sub),    forecastMw, 'k:', 'LineWidth', 2, 'DisplayName', 'Subduction');

axis square; grid on; box on;
legend('Location', 'southeast', 'FontSize', 9);
ylabel('Forecast M_w', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('L_{sub} (km)', 'FontSize', 11, 'FontWeight', 'bold');
ylim([4.5 7.0]);
title('c) Subsection Length (L_{sub})', 'FontSize', 12);
