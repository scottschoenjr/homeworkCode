% Problem 3.7 in Atkinson et al.

clear all
close all
clc

% Set initial conditions
theta0 = [15, 30, 45, 60].*(pi./180);
thetaDot0 = 0;

% Set gravity and arm length
g = 9.81; % [m/s^2]
l = 1.*(0.3048); % [m]

% Set time scale and step
t0 = 0;
tEnd = 10;
dt = 1E-4;

for thetaCount = 1 : length( theta0 )
    
    % Get that initial condition
    y0 = [theta0(thetaCount); thetaDot0];
    
    % Solve system
    [tSol, ySol] = euler_forward( t0, y0, tEnd, dt, 'simplePendulum' );
    
    % Plot result
    theta = ySol(:, 1);
    thetaDot = ySol(:, 2);
    
    figure(1);
    
    subplot( 2, 1, 1 );
    hold all;
    if thetaCount == 1
        set( gca, 'FontSize', 22 );
        ylabel( '$\theta$ [deg]', ...
            'FontSize', 24, 'Interpreter', 'LaTeX' );
        ylim([-75, 75]);
        set( gca, 'YTick', -60:30:60, 'XTickLabel', '' );
        grid on;
    end
    plot( tSol, 180.*theta./pi );
    
    subplot( 2, 1, 2 );
    hold all;
    if thetaCount == 1
        set( gca, 'FontSize', 22 );
        ylabel( '$\dot{\theta}$ [deg/s]', ...
            'FontSize', 24, 'Interpreter', 'LaTeX' );
        xlabel( '$t$ [s]', 'FontSize', 24, 'Interpreter', 'LaTeX' );
        ylim(360.*[-1, 1]);
        set( gca, 'YTick', -360:180:360 );
        grid on;
    end
    plot( tSol, 180.*thetaDot./pi );
    
    % Consider frequency domain signal
    thetaPadded = padarray( theta, 1E6, 0, 'both' );
    dt = tSol(2) - tSol(1);
    Fs = 1./dt;
    fVec = linspace( 0, Fs, length(thetaPadded) ); 
    
    
    figure(2)
    hold all;
    
    if thetaCount == 1
        set( gca, 'FontSize', 22, 'XTick', 0:0.25:5, 'YTickLabel', '' );
        xlabel( '$f/f_{0}$', 'Interpreter', 'LaTeX', 'FontSize', 28 );
        xlim([0.7, 1.3]);
        ylabel( 'Normalized Spectrum' );
        grid on;
    end
    
    fNorm = sqrt(g./l)./(2.*pi);
    s = abs(fft(thetaPadded));
    plot( fVec./fNorm, s./max(s) + 1.0.*(thetaCount - 1) );

    
end



% Define equation of motion
function eom = simplePendulum( odeVars, timeVector, ...
    gravity, pendulumLength )

% Get initial conditions
theta = odeVars(1);
thetaDot = odeVars(2);

y1 = thetaDot;
y2 = -(gravity./pendulumLength).*sin(theta);
eom = [ y1, y2 ];

end



% Forward Euler's method subroutine
function [tVec, y] = euler_forward(t0, y0, tEnd, h, fcn )

% Get order of problem (i.e., number of variables to solve for)
order = length(y0);

% Get number of time points
Nt = floor( (tEnd-t0)./h ) + 1;

% Create time vector for given time step
tVec = linspace(t0, t0 + (Nt-1).*h, Nt)';

% Initialize solution
y = zeros(Nt, order);

% Initial conditions
y(1,:) = y0;

% Get parameter values from workspace
gValue = evalin( 'base', 'g' );
lValue = evalin( 'base', 'l' );

% Compute at subsequent time steps
for tCount = 2:Nt
  y(tCount,:) = ...
      y(tCount-1, :) + h.*feval( fcn, y(tCount-1,:), tVec(tCount-1), ...
      gValue, lValue );
end

end