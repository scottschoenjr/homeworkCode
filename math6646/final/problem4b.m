% Math 6646 - Final Exam
%   PRoble 4b - Shooting BVP

clear all
close all
clc

% Define end points and function values
t0 = 0;
t1 = 1;

x0 = 0;
x1 = 1;

% Set initial guesses for slope
s1 = 0.5;
s2 = 1.5;

% Set number of steps
N = 100;
h = (t1 - t0)./N;

% Solve with the guessed initial conditions
ic1 = [x0; s1];
[t, v] = eulerForward( t0, ic1, t1, h, 'xPrime' );

ic2 = [x0; s2];
[~, w] = eulerForward( t0, ic2, t1, h, 'xPrime' );

% Compute theta [Eqs. (4.3) and (4.2)]
theta = (x1 - v(end, 1))./(w(end, 1) - v(end, 1));
x = theta.*w + (1 - theta).*v;

% Plot
figure()
hold all;

plot( t, x(:, 1), 'k', 'LineWidth', 2 );
plot( t, v(:, 1), ':k', 'LineWidth', 2 );
plot( t, w(:, 1), '--k', 'LineWidth', 2 );

xlabel( '$t$', 'Interpreter', 'LaTeX', 'FontSize', 22 );

lh = legend( '~$x(t)$','~$v(t)$','~$w(t)$', ...
    'Interpreter', 'LaTeX' );
lh.EdgeColor = 'none';

% Format plot
set( gca, 'FontSize', 18, 'FontName', 'Garamond' );

% ========= SUBROUTINES ===========

% Define ODE
function eom = xPrime( odeVars, tVec )

% Get initial conditions
x1 = odeVars(1);
x2 = odeVars(2);

x1Prime = x2;
x2Prime = -x1;
eom = [ x1Prime, x2Prime ];

end

% Forward Euler subroutine
function [tVec, x] = eulerForward(t0, x0, tEnd, h, fcn )

% Get order of problem (i.e., number of variables to solve for)
order = length(x0);

% Get number of time points
Nt = floor( (tEnd-t0)./h ) + 1;

% Create time vector for given time step
tVec = linspace(t0, t0 + (Nt-1).*h, Nt)';

% Initialize solution
x = zeros(Nt, order);

% Initial conditions
x(1,:) = x0;

% Compute at subsequent time steps
for tCount = 2:Nt
    x(tCount,:) = ...
        x(tCount-1, :) + h.*feval( fcn, x(tCount-1,:), tVec );
end

end