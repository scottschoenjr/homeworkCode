% Script to plot basis functions for FEM

clear all
close all
clc

N = 10; % Total nodes
nodesToPlot = [0, 2, 5, 10]; % Which to plot

nodeLables = repmat( {''}, [1, N + 1] );

nValue = 5;
NValue = 10;

% Create x-vector
x = linspace( 0, N, 1E4 );

figure();
hold all;

% Plot each
for nCount = 1 : length( nodesToPlot )
    
    % Get local x
    xLocal = x - nodesToPlot( nCount );
    
    % Define piecewise linear function
    phiN = 0.*x + ...
        ( xLocal > -1 & xLocal <= 0 ).*( ...
           xLocal + 1 ) ...
         + ...
        ( xLocal > 0 & xLocal <= 1 ).*( ...
           1 - xLocal );
       
    % Plot
    plot( x, phiN, 'k' );
    
    % Set label to phi number
    nodeLabels{ nodesToPlot(nCount) + 1 } = ...
        num2str( nodesToPlot( nCount ) );
    
end

% Change the labels for the generaic (nth) and last (Nth) node
nodeLabels{ nValue + 1 } = '$n$';
nodeLabels{ NValue + 1 } = '$N$';

set( gcf, 'Position', [400, 200, 1000, 550] );
set( gca, ...
    'Position', [0.12, 0.2, 0.8, 0.7], ...
    'FontSize', 18, 'XTick', 0:N, 'XTickLabel', nodeLabels );

xlabel( '$x$', 'FontSize', 32 );
ylabel( '$\phi_{n}(x)$', 'FontSize', 28  );

