% Euler's Method

clear all
close all
clc

% Define parameter
y0 = 0;
tMin = 0;
tMax = 10;
hVec = [0.2, 0.1, 0.05];

% Define function and its derivative
yPrime = @(t) cos( t ).^(2);

% Define true solution
tTrue = linspace( tMin, tMax, 1000 );
yTrue = atan(tTrue);

% Perform Euler's method integration

for hCount = 1 : length( hVec )
    h = hVec(hCount);
    
    % Set step and time vector
    tVec = 0 : h : tMax;
    
    % Initialize
    yVec = 0.*tVec;
    yn = y0;
    
    % Compute at each time step
    for tCount = 2 : length(tVec)
        
        yNew = yn + yPrime(yn).*h;
        
        yVec(tCount) = yNew;
        yn = yNew;
        
        
        % Get the error at this point
        [~, tIndex] = min( abs( tTrue - tVec(tCount) ) );
        eVec(tCount) = abs((yTrue(tIndex) - yn));
        
        
    end
    
    % Plot for this step
    figure(1);
    hold all
    plot( tVec, yVec, '-ok', 'MarkerFaceColor', 'k' );
    
    figure(2);
    hold all
    plot( tVec, 100.*eVec, '-ok', 'MarkerFaceColor', 'k' );
    
    if hCount == 1
    legStrings{hCount} = sprintf( ' $h = $%1.2f', h );
    else
        legStrings{hCount} = sprintf( ' %1.2f', h );
    end
    
    % Store to compare
    allResults(hCount).t = tVec;
    allResults(hCount).e = eVec;
   
    
end


figure(1);
plot( tTrue, yTrue, '--k' );
xlabel( '$t$', 'Interpreter', 'Latex', 'FontSize', 32 );
ylabel( '$Y(t)$', 'Interpreter', 'Latex', 'FontSize', 32 );

figure(2);
xlabel( '$t$', 'Interpreter', 'Latex', 'FontSize', 32 );
ylabel( '$100\cdot|\epsilon|$', 'Interpreter', 'Latex', 'FontSize', 32 );
set( gca, 'YTick', 1:10, 'XTick', 0:2:10);
lh = legend( legStrings );
lh.FontSize = 22;
lh.EdgeColor = 'none';