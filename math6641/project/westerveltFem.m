% *************************************************************************
% FEM Solution of Westervelt Equation
%
%   Based on code and work by Blas Dirske "FEM Applied to the 1D Westervelt
%   Equation" (2014). Per TU Delft, assuming CC-BY-SA License.
%   Modifications and content by Scott Schoen Jr under MIT License.
%
% *************************************************************************

clear all; 
close all; 
clc; 

tic;

% Material parameters
c0 = 1500;   % Speed of sound [m/s]
rho0 = 1000; % Density [kg/m^3]
f0 = 1E5;    % Center frequency [Hz]
P0 = 500E3;  % Pressure amplitude [Pa]
beta = 20;   % Parameter of nonlinearity

% Simulation parameters
tolerance = 1E-9; % Iteration tolerance

% Shock formation distance
omega0 = 2.*pi.*f0;
xsh = ( rho0.*c0.^(3) )./( beta.*P0.*omega0 );

% Create 1D mesh
L = xsh + 16.*c0./f0;
n = round(72.*f0.*L/c0);
xVec = linspace(0, L, n+1)';

% Create time vector
dt = (L./n)/(c0.*4);
T = L./c0;
k = round(T./dt)+1;
t = 0:dt:(k-1).*dt; % [s]

% Compute M, K, and C matrices for FEM
[M, S, gM, gS] = ...
    getLinearMatrices(xVec, n, c0);
[C1, C2, C3, gC1, gC2, gC3] = ...
    getNonlinearMatrices( xVec, n, rho0, c0, beta );

% Define source term and its 2nd derivative
[ u0t, u0t_tt ] = pulseWaveform(t, P0, f0);

% Create RHS vector
g = sparse(ones(1,k),1:k,-u0t_tt * gM - u0t * gS,n-1,k,k);

% Initial condition
u1 = zeros(n-1,1);
ut1 = zeros(n-1, 1);

% Solving next time step with first order derivative
u2 = u1 + dt.*ut1;

% ------- Solve Linear Problem -------
tic;

% Intialize
W = zeros(n-1,k);
W(:,1) = u1;
W(:,2) = u2;
W(:,3) = (dt^2 * S + M) \ (2 * M * W(:,2) - M * W(:,1) + dt^2 * g(:,3) );

% Explicit backward-difference
for i = 4:k
    W(:,i) = ( dt^(2)*S + 2*M) \ ...
        ( 5*M*W(:,i-1) - 4*M*W(:, i-2) + M*W(:, i-3) + dt^(2)*g(:,i));
end

% ------- Solve Linear Problem -------

U = zeros(n-1,k);
U(:,1) = u1;
U(:,2) = u2;

% Calculate the third time point using a first order scheme.
for i = 3
    
    z = ( dt^(2)*S + M ) \ ( 2*M*U(:, i-1) - M*U(:,i-2) + dt^(2)*g(:,3) );
    zt = 3*P0*ones(size(z));
    
    while( max( abs(z-zt) ) > tolerance.*P0 )
        
        zt = z;
        
        % Define matrices as in Eq. (4.21)        
        AA = sparse(1:n-1,1:n-1,[u0t(i); zt(1:n-2)],n-1,n-1,n-1);
        AB = sparse(1:n-1,1:n-1,zt,n-1,n-1,n-1);
        AC = sparse(1:n-1,1:n-1,[zt(2:n-1); 0],n-1,n-1,n-1);
        N2 = AA * C1 + AB * C2 + AC * C3; 
        clear AA AB AC;
        
        BA = sparse(1:n-1,1:n-1,([u0t(i); zt(1:n-2)] - [u0t(i-1); ...
            U(1:n-2,i-1)])/dt,n-1,n-1,n-1);
        BB = sparse(1:n-1,1:n-1,(zt - U(:,i-1))/dt,n-1,n-1,n-1);
        BC = sparse(1:n-1,1:n-1,([zt(2:n-1); 0] - [U(2:n-1,i-1); 0])/dt,n-1,n-1,n-1);
        N1 = BA * C1 + BB * C2 + BC * C3;
        clear BA BB BC;
        
        % Update solution
        z = ( dt^(2)*S + M - dt*N1 - N2) \ ...
            ( ...
              ( 2*M - dt*N1 - 2*N2 )*U(:,i-1) + ...
              ( -M + N2 )*U(:,i-2) + ...
              dt^(2)*g(:,i) ...
            );
        
    end
    
    % Store and move on
    U(:,i) = z;
    
end

% Solve all the time points with a second order scheme
for i = 4:k
    z = (dt^2 * S + 2 * M) \ (5 * M * U(:,i-1) - 4 * M * U(:,i-2) + M * U(:,i-3) + dt^2 * g(:,i));
    zt = 3 * P0 * ones(size(z));
    
    % Iteratively solve
    while ( max( abs(z-zt) ) > tolerance.*P0 )
        
        zt = z;
        
        AA = sparse( ...
            1:n-1, 1:n-1, [u0t(i); zt(1:n-2)], n-1, n-1, n-1);
        AB = sparse( ...
            1:n-1, 1:n-1, zt, n-1, n-1, n-1);
        AC = sparse( ...
            1:n-1, 1:n-1, [zt(2:n-1); 0], n-1, n-1, n-1);
        N2 = AA * C1 + AB * C2 + AC * C3; 
        clear AA AB AC;
        
        BA = sparse( ...
            1:n-1, 1:n-1, ...
              ( 3/2 *[ u0t(i); zt(1:n-2) ] - ...
                2*[ u0t(i-1); U(1:n-2,i-1) ] + ...
                1/2*[ u0t(i-2); U(1:n-2,i-2) ] ...
              )/dt, ...
            n-1, n-1, n-1);
        BB = sparse( ...
            1:n-1, 1:n-1, ( 3/2*zt - 2*U(:,i-1) + 1/2*U(:,i-2) )/dt, ...
            n-1, n-1, n-1);
        BC = sparse( ...
            1:n-1, 1:n-1, ...
            ( ...
              3/2*[ zt(2:n-1); 0 ] - ...
              2*[ U(2:n-1,i-1); 0 ] + ...
              1/2*[ U(2:n-1,i-2); 0] ...
            )/dt, ...
            n-1, n-1, n-1);
        N1 = BA * C1 + BB * C2 + BC * C3; 
        clear BA BB BC;
        
        z = ( dt^(2)*S + 2*M - (3/2)*dt*N1 - 2*N2) \ ...
            ( ...
              ( 5*M - 2*dt*N1 - 5*N2)*U(:,i-1) + ...
              (-4*M + (1/2)*dt*N1 + 4*N2)*U(:,i-2) + ...
              ( M - N2)*U(:,i-3) + dt^(2)*g(:,i) ...
             );
    end
    
    U(:,i) = z;
end

% Store FEM solutions
U = [u0t; U; zeros(1,k)]; % Nonlinear solution
W = [u0t; W; zeros(1,k)]; % Linear solution

toc;

%% Plotting

close all;

figure()
hold all;
n = 0;

plotEvery = 2800;
tickFactor = 2;
nT = 1; % Tick inex

for pCount = 1.*plotEvery : plotEvery : length( U(1, :) )

    % Get waveform as a function of space
    waveform = U( :, pCount );
    linearWaveform = W( :, pCount );
    
    plot( 1E2.*xVec, waveform + 1.8.*n.*P0, 'k' );
    n = n + 1;
    
    if mod( n, tickFactor ) == 1;
        yTickVector(nT) = 1.8.*(n - 1).*P0;
        yTickLabels{nT} = [ num2str( round(1E6.*t( pCount )) ), ' \mus' ];
        nT = nT + 1;
    end
    
end
set( gca, 'YTick', yTickVector, ...
    'YTickLabel', yTickLabels, ...
    'FontName', 'Garamond', ...
    'TickLabelInterpreter', 'TeX' );

xlabel( '$x$ [cm]', 'Interpreter', 'LaTeX', 'FontSize', 26 );

box off;

% Plot one specific waveform
figure();
hold on;

spaceIndex = 10800;
plot( 1E2.*xVec, W(:, spaceIndex)./1E3, 'Color', [0, 0, 0, 0.25] );
plot( 1E2.*xVec, U(:, spaceIndex)./1E3, 'Color', [0, 0, 0] );

xlabel( '$x$ [cm]', 'Interpreter', 'LaTeX', 'FontSize', 26 );
ylabel( '$P$ [kPa]', 'Interpreter', 'LaTeX', 'FontSize', 26 );

xlim([40, 55]);

% Plot one specific waveform
figure();
hold on;

Fs = 1./dt;
fVector = linspace( 0, Fs, length(t) );


timeIndex = 3000;

s0Tilde = abs( fft( W(timeIndex, :) ) );
s0Tilde = 20.*log10( s0Tilde./max( s0Tilde ) );

sTilde = abs( fft( U(timeIndex, :) ) );
sTilde = 20.*log10( sTilde./max( sTilde ) );

plot( fVector./f0, s0Tilde, 'Color', [0, 0, 0, 0.25] );
plot( fVector./f0, sTilde, 'Color', [0, 0, 0] );

xlabel( '$f/f_{0}$', 'Interpreter', 'LaTeX', 'FontSize', 26 );
ylabel( 'Magnitude [dB]', 'FontSize', 22 );

xlim([0, 5]);
set( gca, 'XTick', 1:5 );

% Plot level of each harmonic

figure()
s = fft( U );
s0 = max(abs(s(:)));

% Get level of each harmonic
f0Index = find( fVector >= 1.*f0, 1 );
shIndex = find( fVector >= 2.*f0, 1 );
thIndex = find( fVector >= 3.*f0, 1 );
fhIndex = find( fVector >= 4.*f0, 1 );

startIndex = 3000;

f0Level = 20.*log10( abs( s( f0Index, startIndex:end ) ) ./ s0 );
shLevel = 20.*log10( abs( s( shIndex, startIndex:end ) ) ./ s0 );
thLevel = 20.*log10( abs( s( thIndex, startIndex:end ) ) ./ s0 );
fhLevel = 20.*log10( abs( s( fhIndex, startIndex:end ) ) ./ s0 );

hold all;
plot( 1E6.*t(startIndex : end), f0Level, 'k' );
plot( 1E6.*t(startIndex : end), shLevel, 'k' );
plot( 1E6.*t(startIndex : end), thLevel, 'k' );
plot( 1E6.*t(startIndex : end), fhLevel, 'k' );

xlabel( 'Time [$\mu$s]', 'FontSize', 22);
ylabel( 'Normalized Level [dB]', 'FontSize', 22);

xlim([100, 500]);
ylim([-78, 6]);
set( gca, 'YTick', [-75:15:0] );

box off;