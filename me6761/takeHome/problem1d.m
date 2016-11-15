%**************************************************************************
%
% ME-6761 Final Project, Problem 1d
%
%   Plotting mode shapes for 2D waveguide.
%
%     Scott Schoen Jr 20161114
%
%**************************************************************************

clear all
close all
clc

% Channel paramters
H = 30; % [m]
f = 100; % [Hz]
c = 1500; % [m/s]

% Compute first few mode shapes
numModes = 4;
k = 2.*pi.*f./c;

% Create figure
figure()
hold all;
box on;
buffer = 4;
xlim( [ buffer./2, (numModes + 0.5).*buffer ] );
xTicks = zeros(1, numModes);

% Plot each mode
for modeCount = 1:numModes

    % Compute wavenumber
    ky = (modeCount - 0.5).*pi./H;
    
    % Create depth vector
    y = linspace( 0, H, 1000 ); % [m]
    
    % Create the mode shape
    psi = sin( ky.*y );
    
    % Plot
    plot( psi + buffer.*modeCount, y, 'k' );
    
    % Create vector to hold tick positions
    xTicks( modeCount ) = buffer.*modeCount;
    xTicksLabels{ modeCount } = num2str( modeCount );
    
end

% Format plot
xlabel('Mode Number', 'FontSize', 18);
ylabel('Depth [m]', 'FontSize', 18);
set( gca, ...
    'YDir', 'Reverse', ...
    'XTick', xTicks, ...
    'XTickLabel', xTicksLabels ...
    );