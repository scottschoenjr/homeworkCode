% Math 6646 - Final Exam
%   Problem 2b - Predictor-Corrector

clear all
close all
clc

% Define IVP
f = @(t, x) t.^(2);
x0 = 1;

% Set true solution
xTrue = @(t) (1./3).*t.^(3) + 1;

% Set time step and range
tMax = 2;
h = 0.1;
t = ( 0 : h : tMax );

% Initialize
x = 0.*t;
x(1) = x0;

% Compute first four terms from explicit scheme
x(1:4) = xTrue(t(1:4));

% Perform predictor-corrector method
for tCount = 4 : length(t) - 1
    
    n = tCount;
    
    % Get predicted value from explicit scheme
    xPred = ab2(x, t, n, h);
    
    % Use for corrector
    xCorrected = am4(x, t, n, h, xPred);
    
    % Store and loop
    x(tCount + 1) = xCorrected;
    
end

% Plot
figure(1);
hold all;

plot( t, x, 'ko' );
plot( t, xTrue(t), 'k', 'LineWidth', 2 );

% Format plot
set( gca, 'FontSize', 18, 'FontName', 'Garamond' );
ylim([0, 4])

xlabel( '$t$', 'Interpreter', 'LaTeX', 'FontSize', 24  );
ylabel( '$x(t)$', 'Interpreter', 'LaTeX', 'FontSize', 24 );

% 2 step Adams-Bashforth
function x_np1 = ab2(x, t, n, h)

% Define RHS
f = @(t, x) t.^(2);

% Explicit calculation
x_np1 = x(n) + (h./2).*( ...
    3.*f(t(n), x(n)) - ...
    f(t(n-1), x(n-1)) ...
    );

end

% 4 step Adams-Moulton
function x_np1 = am4(x, t, n, h, xpred )

% Define RHS
f = @(t, x) t.^(2);

% Explicit calculation with predicted value
x_np1 = x(n) + (h./720).*( ...
    251.*f( t(n) + h, xpred ) + ...
    646.*f( t(n), x(n) ) - ...
    264.*f( t(n-1), x(n-1) ) + ...
    106.*f( t(n-2), x(n-2) ) - ...
    19.*f( t(n-3), x(n-3) ) ...
    );


end