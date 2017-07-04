%% Kuramoto-Sivashinsky in Fourier Space    
% Domain in real space: [0, L]; computational domain: [-pi, pi].
% Equation in real space for the computational domain:
% u_t = -s^4*u_{xxxx} - s^2*u_{xx} - s*stab*(u.^2)_x/2,
% where s = 2*pi/L is the spatial scaling parameter and stab is 
% the energy transfer parameter. 
%
% Equation in Fourier space:
% u_hat_t = symb.*u_hat - q.*( u_hat^{*2} ),
% where symb = s^2*k.^2 - s^4.*k.^4 is the symbol of the linear part of 
% the pde, q = 1i*s*k*stab/2 is the coefficient of the nonlinear part and
% u_hat^{*2} denotes the second convolution power, i.e. convolution of 
% u_hat with itself. 
%
% Temporal discretization: fourth order exponential Runge-Kutta.
%% Physical Parameters                      
L = 13;             % Domain size for the equation before rescaling  
Tfinal = 5e2;       % Total length of the simulation
s = 2*pi/L;         % Scaling parameter
%% Computational Parameters                 
MakeMovie = 0;      % Controls creation of .avi file 
decayRates = 0;     % Controls computation of decay rates of the Fourier coefficients
SN = 64;            % Number of grid points (number of computed modes)
M = 64;             % Number of points for complex means
dx = 2*pi/SN;       % Spatial resolution
x = -pi:dx:pi-dx;   % Physical space
TN = 5e4;           % Number of time steps
plotgap = 1e1;      % Number of time steps between plots
dt = Tfinal/TN;     % Size of the time step 
%% Initial Condition                        
u = sin(2*x);        % Initial condition in real space                   
u_hat = fft(u);      % Initial condition in Fourier space
%% Auxiliary Variables                      
numplots = TN/plotgap;           % Number of plots
data = zeros(numplots+1, SN);    % Preallocating for solution
data(1, :) = u_hat;              % Storing initial condition
Mds = [0 -SN/2+1:SN/2-1];        % Modes ordered in a reasonable way 
Mds = ifftshift(Mds);            % Modes ordered in the Matlab way
K = Mds*s;                       % Wave numbers
evals = K.^2-K.^4;               % Eigenvalues of the linear part of the pde
E1 = exp(dt*evals);              % Full linear step
E2 = exp(dt*evals/2);            % Half linear step
q = -1i*K/2;                     % Coefficient of the nonlinear part of the pde
r = exp(1i*pi*((1:M)-.5)/M);     % Roots of unity
%% Time Stepping Coefficients               
CC = dt*(Mds(:,ones(M,1)) + r(ones(SN,1),:))';  % Complex countours for reciprocal evaluation
Q  = dt*real(mean((exp(CC/2)-1)./CC)); 
f1 = dt*real(mean((-4-CC+exp(CC).*(4-3*CC+CC.^2))./CC.^3)); 
f2 = dt*real(mean((2+CC+exp(CC).*(-2+CC))./CC.^3)); 
f3 = dt*real(mean((-4-3*CC-CC.^2+exp(CC).*(4-CC))./CC.^3));  
%% Time Integration                         
tic
for pic = 2:numplots+1            % Stepping from one plot to the next
    for step = 1:plotgap          % Stepping between the plots
        % Fourth order exponential Runge-Kutta 
        N1 = q.*fft(ifft(u_hat,'symmetric').^2);    % Telling ifft that my vector
        A = E2.*u_hat + Q.*N1;                      % is the Fourier transform of 
        N2 = q.*fft(ifft(A,'symmetric').^2);        % a real function
        B = E2.*u_hat + Q.*N2;
        N3 = q.*fft(ifft(B,'symmetric').^2);
        C = E2.*A + Q.*(2*N3-N1);
        N4 = q.*fft(ifft(C,'symmetric').^2);
        u_hat = E1.*u_hat + N1.*f1 + 2*(N2+N3).*f2 + N4.*f3;
    end
    data(pic, :) = u_hat;    % Storing the solution
end
fprintf('Done with time itegration after %0.0f seconds, now preparing plots.\n', toc);
%% Post-processing of Data                  
Rdata = ifft(data, [], 2, 'symmetric');            % Going back to the real spaces
Du = ifft(data*diag(1i*K), [], 2, 'symmetric');    % Derivative of a solution
pSpec = abs(data(:, 1:SN/2).^2);                   % Power spectrum
means = mean(Rdata, 2);                            % Means at each time step
maxs = max(abs(Rdata), [], 2);                     % Maximums at each time step
L2norms = sqrt(sum(Rdata.^2, 2)*L/SN);             % L2 norm at each time step
time = 0:dt*plotgap:Tfinal;                        % Creating the time vector
%% Plotting                                 
if exist('../KS_Pictures_Movies','dir')~=7
    mkdir ../KS_Pictures_Movies    % Creates folder for outputs
    disp('Expected folder for outputs was not found, so I made it.')
end  
fig100 = figure(100);
set(fig100, 'PaperOrientation', 'landscape');
set(fig100, 'position', [0 0 1280 800]);
p1 = subplot(2, 2, 1);
pos1 = get(p1, 'position');
pos1 = pos1 + [-.05 -.1 .1 .1]; 
set(p1, 'position', pos1);
surf((x+pi)/s, time, Rdata, 'edgecolor', 'none')   
colorbar, view([0 90])
title(['Time evolution in real space with ETDRK4, dt = ', num2str(dt)], 'fontsize', 16)
xlim([0 L-dx/s]); ylim([0 Tfinal]);
xlabel('Space', 'Fontsize', 16), ylabel('Time', 'Fontsize', 16)
colormap(p1, 'default') 

p2 = subplot(2, 2, 2);
pos2 = get(p2, 'position');
pos2 = pos2 + [-.02 -.1 .1 .1]; 
set(p2, 'position', pos2);
contourf(K(1:SN/8), time, log(pSpec(:, 1:SN/8)), 30, 'edgecolor', 'none')
colormap(p2, autumn), colorbar
title('Power spectrum evolution on logarithmic scale', 'fontsize', 16) 
xlabel('Wave numbers', 'Fontsize', 16), ylabel('Time', 'Fontsize', 16)

p3 = subplot(2, 1, 2);
pos3 = get(p3, 'position');
pos3 = pos3 + [0 0 0 -.1]; 
set(p3, 'position', pos3);
plot(time, L2norms, '.', 'markersize', 15)
hold on
plot(time, maxs, '.', 'markersize', 15), hold off
title('Time evolution of norms', 'fontsize', 16)
xlabel('Time', 'Fontsize', 16)
legend({'L^2 norm of u', 'Maximum of u'}, 'Location', 'Northwest', 'FontSize', 12)

print(fig100, sprintf('../KS_Pictures_Movies/KS_L%0.2f.png', L), '-dpng')

% figure(101)
% surf((x+pi)/s, time, Du, 'edgecolor', 'none')
% colorbar, view([0 90])
% title('Time evolution of derivative')
% xlim([0 L-dx/s]); ylim([0 Tfinal]);
% xlabel('Space', 'Fontsize', 16), ylabel('Time', 'Fontsize', 16)
%% Movie                                    
if MakeMovie
name = sprintf('../KS_Pictures_Movies/KS_modes_L%0.2f.avi', L);  % Name for the video file
vidObj = VideoWriter(name);                                      % Creating video file
vidObj.FrameRate = 10;  open(vidObj);                            % Opening video file
         
fig102 = figure(102);                           % Initial plots
set(fig102, 'position', [0 0 1280 800]);        % Full screen figure
plot102a = plot(K(1:SN/2), log(pSpec(1, 1:SN/2)),...
    '-ro', 'MarkerFaceColor', 'r',...
    'LineWidth', 1, 'MarkerSize', 8); 
hold on
%plot102b = plot(K(1:SN/8), 0*(1:SN/8), 'linewidth', 2);
plot102c = plot(K(1:SN/2), 60*evals(1:SN/2), 'k--',  'linewidth', 1);
hold off

gca102=gca; axis([K(1) K(SN/2) -80 20]), grid on
title('Time = 0', 'Fontsize', 20)
xlabel('Wave numbers', 'Fontsize', 16)
ylabel('Energy', 'Fontsize', 16)
writeVideo(vidObj, getframe(fig102));
tic
for t = 2:numplots+1

    %rt = polyfit(10:SN/8, log(pSpec(t, 10:SN/8)), 1);  % Best fit linear model
    %rates(t) = -rt(1);                                 % Convergence rate
    
    set(plot102a,'YData', log(pSpec(t, 1:SN/2)))
    %set(plot102b,'YData', polyval(rt, 1:SN/8))
    %title(gca102, sprintf('Time = %1.2f, decay rate = %1.2f',...
        %(t-1)*plotgap*dt, -rt(1)), 'Fontsize', 14)
    title(gca102, sprintf('Time = %0.1f', (t-1)*plotgap*dt), 'Fontsize', 20)
    %legend([plot102a, plot102b], {'Power Spectrum', 'Best Fit Linear Model'}, 'Fontsize', 14)
    legend([plot102a, plot102c], {'Power Spectrum', 'Eigenvalues'}, 'Fontsize', 16)
    writeVideo(vidObj, getframe(fig102));
    
end, toc
close(vidObj); 
end
%% Convergence Rate Evolution               
if decayRates 
rates = zeros(1, numplots+1);      % Preallocating for convergence rates
figure(103)
plot(numplots/4:numplots+1, rates(numplots/4:end), 'linewidth', 3)
xlim([numplots/4 numplots+1])
title(['Exponential rates, mean = ', num2str( mean(rates(numplots/4:end)) )], 'Fontsize', 14) 
ww = fft(rates(numplots/4:end));
figure(104)
plot(2:numplots/2, abs(ww(2:numplots/2)), '.', 'markersize', 15)
title('Frequencies of decay rates oscillations')
end
%% Notes                                    
% L=40,   u0 = cos(3x) gives steady state solution for t<100
% L=50,   u0 = 2.5cos(6x) gives steady state solution for t<1e+4
% L=500,  u0 = cos(6x) gives cells
% L=2\pi, u0 = exp(-10*x.^2), antidiffusivity = 6 gives travelling waves