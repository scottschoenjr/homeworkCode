% Problem 4

clear all
close all
clc

% Define parameters
m = 4;
c = 1;
k = 2;

% Define shorthands
zeta = c./(2.*m);
omega0 = sqrt(k./m);

% Define initial conditons
x0 = 1;
v0 = -8;

% Define analytical solution
t = linspace( 0, 20, 1000 );
omega1 = abs( sqrt( zeta^(2) - omega0.^(2) ) );
C1 = 1;
C2 = ( -8 + 1./8 )./omega1;
xAnalytical = exp( -zeta.*t ).*( ...
    C1.*cos( omega1.*t ) + C2.*sin( omega1.*t ) ...
    );

% Plot analytical solution


% Compute numerical solution as a check
eom = @(tVar, xVar)  [xVar(2); -(c./m).*xVar(2) - (k./m).*xVar(1)];
initialConditions = [x0; v0]; % Initial conditions
tRange = [min(t), max(t)];
tSpan = max(t) - min(t);
odeOptions = odeset( ...
    'reltol', 1E-6, ...
    'abstol', 1E-6, ...
    'maxstep', tSpan./100 ...
    );
[ tNum, xNum ] = ode45( eom, tRange, initialConditions, odeOptions);

figure(1)
hold on;
plot( t, xAnalytical, 'k' );
plot( tNum, xNum(:, 1), 'ko', 'MarkerSize', 3);
xlabel( 'Time [s]', 'FontSize', 20 );
ylabel( 'Position [ft]', 'FontSize', 20  );


% Add legend
legH = legend( '~Eq. (4.11)', '~Numerical' );
legH.FontSize = 20;

    

