% Math 6646 - Final Exam
%   Problem 1 - Forward, backward, trapezoidal comparison

clear all
close all
clc

% Set initial condition
y0 = 1;

% Set time scale and step
t0 = 0;
dt = 0.005;
tEnd = 3.5;

% Set tolerance for implicit method
tolerance = 0.001;
maxIterations = 20;

% Solve system with each method

% Forward Euler
[tSol, ySol_fe] = ...
    eulerForward( t0, y0, tEnd, dt, 'yPrime' );

% Backward Euler
[~, ySol_be] = ...
    eulerBackward( t0, y0, tEnd, dt, 'yPrime', tolerance, maxIterations );

% Trapezoidal rule
[~, ySol_tr] = ...
    trapezoidalRule( t0, y0, tEnd, dt, 'yPrime', tolerance, maxIterations );

% Plot result
figure();
hold all;
set( gca, 'FontSize', 18, 'FontName', 'Garamond', 'XTick', 0:5 );

expSolPlot = plot( tSol, ySol_fe, '--k', 'LineWidth', 2 );
impSolPlot = plot( tSol, ySol_be, 'k', 'LineWidth', 2 );
trSolPlot = plot( tSol, ySol_tr, ':k', 'LineWidth', 2 );

% axis equal;
box off;

ylim([-0.75, 1.25]);
xlim([0, 3.5]);

xlabel( '$t$', 'FontSize', 24, 'Interpreter', 'LaTeX' );
ylabel( '$y$', 'FontSize', 24, 'Interpreter', 'LaTeX' );

lh = legend( ' Forward', ' Backward', ' Trapezoidal' );
lh.EdgeColor = 'none';

% ======= SUBROUTINES =======

% Define RHS
function f = yPrime( tVec, y )

% Get initial conditions
f = -100.*(y - sin(tVec));

end

% Forward Euler method
function [tVec, y] = eulerForward(t0, y0, tEnd, h, fcn )

% Get number of time points
Nt = floor( (tEnd-t0)./h ) + 1;

% Create time vector for given time step
tVec = linspace(t0, t0 + (Nt-1).*h, Nt)';

% Initialize solution
y = zeros(Nt, 1);

% Initial condition
y(1) = y0;

% Compute at subsequent time steps
for tCount = 2:Nt
    
    tn = tVec(tCount - 1);
    yn = y(tCount - 1);
    fn = feval( fcn, tn, yn );
    
    y(tCount) = y(tCount-1) + h.*fn;
end

end

% Backward Euler method
function [t,ySol] = eulerBackward(t0, y0, tEnd, h, fcn, ...
    tolerance, maxIterations)

% Set defaults for optional variables
if nargin < 6
    tolerance = 1E-4;
end
if nargin < 7
    maxIterations = 20;
end

% Get number of time points
Nt = floor( (tEnd-t0)./h ) + 1;

% Create time vector for given time step
t = linspace(t0, t0 + (Nt-1).*h, Nt)';

ySol = zeros(Nt, 1);
ySol(1) = y0;
tCount = 2;


while tCount <= Nt

    % Get starting estimate with Forward Euler
    yt1 = ySol(tCount-1, :) + ...
        h.*feval(fcn, t(tCount-1), ySol(tCount-1));
    
    % Fixed point iteration
    
    % Initialize
    iterationCount = 0;
    diff = 2.*tolerance;
    
    % Iterate until convergence
    while ( diff > tolerance ) && ( iterationCount < maxIterations )
        
        % Iterate
        yt2 = ySol(tCount-1, :) + h*feval(fcn, t(tCount), yt1);
        
        % Get max difference
        diff = max(abs( yt2-yt1) );
        
        % Update starting value
        yt1 = yt2;
        
        % Increment iteration
        iterationCount = iterationCount +1;
    end
    
    % Warn user if it did not converge
    if iterationCount >= maxIterations
        warningString = sprintf( ...
            [ 'Fixed-point did not converge below %4.2f ', ...
              'after %3d iterations [residual = %4.2f].'], ...
            tolerance, maxIterations, diff );
        warning( warningString ); 
    end
    
    % Store solution and go to next time step
    ySol(tCount) = yt2;
    
    tCount = tCount+1;
end

end


% Trapezoidal rule
function [t,ySol] = trapezoidalRule(t0, y0, tEnd, h, fcn, ...
    tolerance, maxIterations)


% Set defaults for optional variables
if nargin < 6
    tolerance = 1E-4;
end
if nargin < 7
    maxIterations = 20;
end

% Get number of time points
Nt = floor( (tEnd-t0)./h ) + 1;

% Create time vector for given time step
t = linspace(t0, t0 + (Nt-1).*h, Nt)';

ySol = zeros(Nt, 1);
ySol(1) = y0;
tCount = 2;


while tCount <= Nt
    
    % Get forward Euler slope
    fn = feval(fcn, t(tCount-1), ySol(tCount-1));

    % Get starting estimate with Forward Euler
    yt1 = ySol(tCount-1, :) + h.*fn;
    
    % Fixed point iteration
    
    % Initialize
    iterationCount = 0;
    diff = 2.*tolerance;
    
    % Iterate until convergence
    while ( diff > tolerance ) && ( iterationCount < maxIterations )
        
        % Iterate
        yt2 = ySol(tCount-1, :) + h*feval(fcn, t(tCount), yt1);
        
        % Get max difference
        diff = max(abs( yt2-yt1) );
        
        % Update starting value
        yt1 = yt2;
        
        % Increment iteration
        iterationCount = iterationCount +1;
    end
    
    % Warn user if it did not converge
    if iterationCount >= maxIterations
        warningString = sprintf( ...
            [ 'Fixed-point did not converge below %4.2f ', ...
              'after %3d iterations [residual = %4.2f].'], ...
            tolerance, maxIterations, diff );
        warning( warningString ); 
    end
    
    % Do final iteration to get backward Euler slope
    fnp1 = feval(fcn, t(tCount), yt1);
    
    % Store trapezoidal solution and go to next step
    ySol(tCount) = ySol(tCount-1, :) + (h./2).*( fn + fnp1 );
    
    tCount = tCount+1;
end


end
