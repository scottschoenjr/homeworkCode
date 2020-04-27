% Math 6646 - Final Exam
%   Problem 3b - RK3

clear all
close all
clc

% Set free parameter
alpha = 0.75;

% Define IVP
f = @(t, x) t.^(2);
x0 = 1;

% Set true solution
xTrue = @(t) (1./3).*t.^(3) + 1;

% Set time step and range
tMax = 1;
hVec = [0.2, 0.1, 0.05, 0.025];

% Define Butcher table constants
c1 = 0;
c2 = alpha;
c3 = 0.5;

b1 = 0;
b2 = 0;
b3 = 1;

a21 = alpha;
a31 = 0;
a32 = 0.5;

% Integrate at each step for each step size
for hCount = 1 : length(hVec)
    
    % Assemble time vector for this step size
    h = hVec(hCount);    
    t = ( 0 : h : tMax );
    x = 0.*t;
    x(1) = x0;
    
    
    for count = 2 : length(t)
        
        xn = x(count - 1);
        tn = t(count);
        
        k1 = f(tn, xn);
        k2 = f(tn + c2.*h, xn + h.*a21.*k1);
        k3 = f(tn + c3.*h, xn + h.*a31.*k1 + h.*a32.*k2);
        
        xnp1 = xn + h.*(b1.*k1 + b2.*k2 + b3.*k3);
        
        x(count) = xnp1;
        
    end
    
    % Compute error
    e = x - xTrue(t);
    
    % Plot
    figure(1);
    hold all;
    
    plot( t, abs(e), 'k', 'LineWidth', 2 );
    
    % Save
    allErrors(hCount).e = e;
    allErrors(hCount).t = t;
    
end

% Format plot
set( gca, 'FontSize', 18, 'FontName', 'Garamond', 'YScale', 'log', ...
    'YTick', [1E-6, 1E-4, 1E-2, 1], 'YTickLabel', ...
    {'10^{-6}'; '10^{-4}'; '0.01'; '1'} );

xlabel( '$t$', 'Interpreter', 'LaTeX', 'FontSize', 24  );
ylabel( '$|x(t) - h_{h}(t)|$', 'Interpreter', 'LaTeX', 'FontSize', 24 );