% Math 6646 - Final Exam
%   Finite difference BVP [uses code from HW4]

clear all
close all
clc

% Define (normalized) coefficient functions
p = @(t) 0.*t;
q = @(t) 0.*t - 1;
r = @(t) 0.*t;

% Define end points and function values
t0 = 0;
t1 = 1;

x0 = 0;
x1 = 1;

% Set number of steps
N = 100;
h = (x1 - x0)./N;

% Initialize solution vectors
tSolution = linspace( t0, t1, N + 1)';
xSolution = zeros( N + 1, 1 );

% Enforce BCs
xSolution(1) = x0;
xSolution(end) = x1;

% Evaluate coefficient functions at each time t
pVals = p( tSolution );
qVals = q( tSolution );

% Create tridiagonal matrix
a =  0.5.*h.*pVals(2:N-1) - 1;
b =  2 + h.^(2).*qVals(2:N);
c = -0.5.*h.*pVals(3:N) - 1;
A = diag(a, 1) + diag(b) + diag(c, -1);

% Create RHS vector
rVals = r( tSolution );
f = -h.^(2).*rVals(2:N);
f(1) = f(1) + (1 + (h./2).*pVals(2) ).*x0;
f(end) = f(end) + (1 - (h./2).*pVals(N) ).*x1;

% Solve
xSolution(2:N) = A\f;

% Plot
figure();
hold all;

plot( tSolution, xSolution, 'k', 'LineWidth', 2 );

% Format plot
xlabel( '$t$', 'Interpreter', 'LaTeX', 'FontSize', 22 );
ylabel( '$x(t)$', 'Interpreter', 'LaTeX', 'FontSize', 22 );
set( gca, 'FontSize', 18, 'FontName', 'Garamond' );
