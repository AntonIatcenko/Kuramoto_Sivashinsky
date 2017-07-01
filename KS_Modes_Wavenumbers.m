%% Modes and Wavenumbers Diagramm      
% for the linear part of the 
% Kuramoto-Sivashinsky equation
%% Parameters                          
L = linspace(000, 100);                  % Possible domain size
modes1 = L/(2*sqrt(2)*pi);               % Dominating modes
mmin1 = ceil(modes1(1));                 % Smallest interger mode
mmax1 = floor(modes1(end));              % Largest integer mode
intmodes1 = mmin1:mmax1;                 % Integer dominating modes
L_intmodes1 = 2*sqrt(2)*pi*intmodes1;    % Corresponding domain sizes

modes2 = L/(2*pi);                       % Nontrivial steady modes
mmin2 = ceil(modes2(1));                 % Smallest interger mode
mmax2 = floor(modes2(end));              % Largest integer mode
intmodes2 = mmin2:mmax2;                 % Integer dominating modes
L_intmodes2 = 2*pi*intmodes2;            % Corresponding domain sizes

%% Plot                                
figure(110)
plot(L, modes1, 'b', 'linewidth', 2), grid on, hold on
plot(L_intmodes1, intmodes1, 'r.', 'markersize', 30)
plot(L, modes2, 'm', 'linewidth', 2)
plot(L_intmodes2, intmodes2, 'k.', 'markersize', 30);
hold off
title('Modes and Wavenumbers Diagramm', 'fontsize', 20)
xlabel('Domain Size', 'fontsize', 16)
ylabel('Mode Number', 'fontsize', 16)
legend({'Dominating modes', 'Integer dominating modes',...
    'Nontrivial neutral modes', 'Integer nontrivial neutral modes'},...
    'location', 'northwest', 'fontsize', 16)
set(gca, 'ytick', intmodes2)
%set(gca, 'xtick', round(L_intmodes1, 2))
set(gca, 'xtick', round(L_intmodes2, 2))