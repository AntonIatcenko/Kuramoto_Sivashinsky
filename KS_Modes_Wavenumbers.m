%% Modes and Wavenumbers Diagramm                 
% for the linear part of the 
% Kuramoto-Sivashinsky equation
%% Parameters for Comparison Across Domains       
L = linspace(000, 100, 500);             % Possible domain size
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
%% Plot Across Domains                            
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
set(gca, 'xtick', round(L_intmodes1, 2))
%% Parameters for a Single Domain                 
L = 13;                       % Size of the domain
s = 2*pi/L;                   % Scaling parameter
n = 6;                        % Number of eigenvalues to be plotted
xmax = n*s + 1/4;             % Range for the plot
x = linspace(0, xmax, 500);   % Plotting domain
evals = x.^2 - x.^4;          % Eigenvalues of the linear part of KS
mds = s*(0:n);                % Modes to be shown
mdsevals = mds.^2 - mds.^4;   % Corresponding eigenvalues
%% Plot for a Single Domain                       
figure(111)
plot(x, evals, 'b-', 'linewidth', 2), hold on
plot(mds, mdsevals, 'r.', 'markersize', 40)
plot(x, 0*x, 'k', 'linewidth', 1), hold off
title(['Modes for domain of size ', num2str(L)], 'fontsize', 20)
axis([0 xmax min(evals) max(evals)-min(evals)/4]), grid on
xticks(0:0.2:mds(end)+0.2), xticklabels(string(0:0.2:n))