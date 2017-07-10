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
%% Spatial Parameters                       
L = 13.3;           % Domain size for the equation before rescaling 
s = 2*pi/L;         % Scaling parameter
SN = 512;           % Number of grid points (number of computed modes)
dx = 2*pi/SN;       % Spatial resolution
x = -pi:dx:pi-dx;   % Physical space
%% Temporal Parameters                      
Tfinal = 2e2;                       % Total length of the simulation
plotgap = 1;                        % Physical time between subsequent plots
dt = 1e-3;                          % Initial size of a time step
plTN = ceil(plotgap/dt);            % Number of time steps between plots
dt = plotgap/plTN;                  % Adjusted time step
numplots = ceil(Tfinal/plotgap);    % Number of plots
TN = plTN*numplots;                 % Total number of time steps
%% Initial Conditions                       
u = sin( (1:4)'.*x );         % Initial condition in real space                   
u_hat = fft(u, SN, 2);        % Initial condition in Fourier space
%% Output Parameters                        
ICnum = numel(u(:, 1));       % Number of initial conditions
ICplots = 1:ICnum;            % Solutions to be plotted
%% Preallocations                           
data = zeros(numplots+1, ICnum, SN);    % Preallocating for solution
data(1, :, :) = u_hat;                  % Storing initial condition
%% Fourier Variables                        
Mds = [0 -SN/2+1:SN/2-1];   % Modes ordered in a reasonable way 
Mds = ifftshift(Mds);       % Modes ordered in the Matlab way
K = Mds*s;                  % Wave numbers
evals = K.^2-K.^4;          % Eigenvalues of the linear part of the pde
%% Time Stepping Coefficients               
E1 = exp(dt*evals);   E2 = exp(dt*evals/2);       % Full and half linear steps
q = -1i*K/2;                                      % Coefficient of the nonlinear part of the pde
M = 64;                                           % Number of points for complex means
r = exp(1i*pi*((1:M)-.5)/M);                      % Roots of unity
CC = dt*(Mds(:,ones(M,1)) + r(ones(SN,1),:))';    % Complex countours for reciprocal evaluation
Q  = dt*real(mean((exp(CC/2)-1)./CC)); 
f1 = dt*real(mean((-4-CC+exp(CC).*(4-3*CC+CC.^2))./CC.^3)); 
f2 = dt*real(mean((2+CC+exp(CC).*(-2+CC))./CC.^3)); 
f3 = dt*real(mean((-4-3*CC-CC.^2+exp(CC).*(4-CC))./CC.^3)); 
%% Time Integration                         
tic, fprintf('Starting time integration, going to make %2.0e steps ... ', TN)
waitbarhandle = waitbar(0, 'Starting time integration');
titleHandle = get(findobj(waitbarhandle , 'Type','axes'),'Title');
set(titleHandle, 'FontSize', 24)
for tcount = 1:TN
    % Fourth order exponential Runge-Kutta 
    N1 = q.*fft(ifft(u_hat, SN, 2, 'symmetric').^2, SN, 2);    % Telling ifft that my vector
    A = E2.*u_hat + Q.*N1;                                     % is the Fourier transform of 
    N2 = q.*fft(ifft(A, SN, 2, 'symmetric').^2, SN, 2);        % a real function
    B = E2.*u_hat + Q.*N2;
    N3 = q.*fft(ifft(B, SN, 2, 'symmetric').^2, SN, 2);
    C = E2.*A + Q.*(2*N3-N1);
    N4 = q.*fft(ifft(C, SN, 2, 'symmetric').^2, SN, 2);
    u_hat = E1.*u_hat + N1.*f1 + 2*(N2+N3).*f2 + N4.*f3;
    if mod(tcount, plTN) == 0
    data(tcount/plTN+1, :, :) = u_hat;    % Storing the solution
    waitbar(tcount/TN, waitbarhandle,...
        sprintf('Time %0.0f out of %0.0f', dt*tcount, Tfinal))
    end
end
close(waitbarhandle)
fprintf('done after %0.0f seconds, now preparing plots.\n', toc)
%% Post-processing the Data                 
Rdata = ifft(data, SN, 3, 'symmetric');            % Going back to the real spaces
pSpec = abs(data(:, :,  1:SN/2).^2);               % Power spectrum
means = mean(Rdata, 3);                            % Means at each time step
maxs = max(abs(Rdata), [], 3);                     % Maximums at each time step
L2norms = sqrt(sum(Rdata.^2, 3)*L/SN);             % L2 norm at each time step
time = 0:plotgap:Tfinal;                           % Creating the time vector
%% Plotting                                 
fig0 = figure(500);
set(fig0, 'PaperOrientation', 'landscape');
set(fig0, 'position', [0 0 1280 800]);
for runNum = 1:ICnum
j = ICplots(runNum);
fig = figure(500 + j);
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'position', [0 0 1280 800]);
surf((x+pi)/s, time, squeeze(Rdata(:, j, :)), 'edgecolor', 'none')   
colorbar, view([0 90])
title(['Time evolution for initial condition ', num2str(j)], 'fontsize', 20)
xlim([0 L-dx/s]); ylim([0 Tfinal]);
xlabel('Space', 'Fontsize', 16), ylabel('Time', 'Fontsize', 16) 

figure(500)
spl = subplot(ICnum, 1, j);
plot(time, L2norms(:, j), '.', 'markersize', 15), hold on
plot(time, maxs(:, j), '.', 'markersize', 15), hold off
title(['Time evolution of norms for initial condition ', num2str(j)], 'fontsize', 20)
xlabel('Time', 'Fontsize', 16)
legend({'L^2 norm of u', 'Maximum of u'}, 'Location', 'Northwest', 'FontSize', 12)

end