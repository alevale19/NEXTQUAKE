forecastMw = 5.5:0.1:8;

% Inizializziamo i vettori e le matrici per ospitare i risultati del ciclo
numSteps = length(forecastMw);
av = zeros(1, numSteps);

% Pre-allocazione per le matrici dei modelli rimanenti (solo Rupture Area)
% Gruppo WC94 RA: contiene 4 relazioni (N, SS, R)
L_WC94_all = zeros(3, numSteps);
% Gruppo T17 RA Crustal: contiene 3 relazioni (N, SS, R) -> Esclusa la subduzione
L_T17_all  = zeros(3, numSteps);

% Generiamo i dati ciclando magnitudo per magnitudo
for i = 1:numSteps
    % Utilizziamo 'evalc' per silenziare i log della funzione durante il ciclo for
    % in modo da non intasare la riga di comando di MATLAB
    [txt, subSectLength] = evalc("seismicScalingEngine(forecastMw(i), 'Mw2Length', 'All')");
    av(i) = subSectLength;
    
    % Gruppo 1: Wells & Coppersmith (1994) Rupture Area -> sqrt(RA)
    L_WC94_all(:, i) = [sqrt(10^(-2.87 + 0.82 * forecastMw(i))); ...
                        sqrt(10^(-3.42 + 0.90 * forecastMw(i))); ...
                        sqrt(10^(-3.99 + 0.98 * forecastMw(i)))];
                    
    % Gruppo 2: Thingbaijam (2017) Crustal Rupture Area -> sqrt(RA)
    L_T17_all(:, i)  = [sqrt(10^(-2.551 + 0.808 * forecastMw(i))); ...
                        sqrt(10^(-3.486 + 0.942 * forecastMw(i))); ...
                        sqrt(10^(-4.362 + 1.049 * forecastMw(i)))];
end

% =====================================================================
% GENERAZIONE DEL GRAFICO INFORMATIVO (Solo Modelli di Area)
% =====================================================================
figure('Name', 'Scaling Relations Calibration (RA Only)', 'Color', 'w');

% Plot delle singole relazioni (passate trasposte per plottare le curve)
semilogy(forecastMw, L_WC94_all', 'r', 'LineWidth', 1)
hold on;
semilogy(forecastMw, L_T17_all', 'b', 'LineWidth', 1)

% Plot della nuova media finale (l'ensemble ridotto "av") in forte evidenza
semilogy(forecastMw, av, 'k', 'LineWidth', 3) 

% Linea di controllo a 15 km (limite classico UCERF3/Appennino per le sezioni)
plot([5.5, 8.0], [15, 15], 'g--', 'LineWidth', 2) 

xlabel('Magnitude (Mw)');
ylabel('Length / Sqrt(Area) (km)');
title('Ensemble Scaling Relations for Subsection Sizing (Area Models Only)');

% Sistemazione della legenda per l'assetto aggiornato (7 curve totali + 2 evidenze)
legend('WC94 RA Models', '', '', '', ...
       'T17 RA Crustal Models', '', '', ...
       'FINAL ENSEMBLE MEAN', 'Target Limit (15 km)', 'Location', 'Best');
   
grid on;