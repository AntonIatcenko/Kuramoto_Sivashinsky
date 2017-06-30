%% Modes and Wavenumbers Diagramm      
% for the linear part of the 
% Kuramoto-Sivashinsky equation
%% Parameters                          
L = linspace(000, 100);                 % Possible domain size
modes = L/(2*sqrt(2)*pi);               % Dominating modes
mmin = ceil(modes(1));                  % Smallest interger mode
mmax = floor(modes(end));               % Largest integer mode
intmodes = mmin:mmax;                   % Integer dominating modes
L_intmodes = 2*sqrt(2)*pi*intmodes;     % Corresponding domain sizes
%% Plot                                
figure(110)
plot(L, modes, 'linewidth', 2), grid on, hold on
pint = plot(L_intmodes, intmodes, 'r.', 'markersize', 30);
hold off
title('Modes and Wavenumbers Diagramm', 'fontsize', 20)
xlabel('Domain Size', 'fontsize', 16)
ylabel('Dominant Mode Number', 'fontsize', 16)
legend({'All modes', 'Integer modes'},...
    'location', 'northwest', 'fontsize', 16)
set(gca, 'xtick', round(L_intmodes, 2)), set(gca, 'ytick', intmodes)