%**************************************************************************
%
% ME-6761 Final Project, Problem 4c
%
%   Computing phase speed plot
%
%     Scott Schoen Jr 201611130
%
%**************************************************************************

clear all
close all
clc

% Channel paramters
H = 30; % [m]
fVector = linspace(1, 500, 1E5); % [Hz]
c = 1500; % [m/s]

% Compute first few mode shapes
numModes = 4;

% Initialize
c_gr = zeros( 1, length(fVector) );
c_ph = zeros( 1, length(fVector) );

% Create figure
figure()
set( gcf, 'Position', [100, 100, 1200, 800], 'Color', 'w' );
set( gca, 'FontSize', 18 );
box on;
hold all

% Plot each mode
for modeCount = 1:numModes
    
    % At each frequency
    for fCount = 1:length(fVector)
        
        % Get wavenumber
        omega = 2.*pi.*fVector(fCount);
        k = omega./c;
        
        % Compute wavenumber
        ky = (modeCount - 0.5).*pi./H;
        
        % Compute phase speed
        c_gr(fCount) = c.*(1 - (ky./k).^(2)).^(0.5);
        
        % Compute phase speed
        c_ph(fCount) = c.^(2)./c_gr(fCount);
        
    end
    
    % Plot NaNs for imaginary parts
    nanInds = find( imag( c_ph ) ~= 0 );
    c_ph( nanInds ) = NaN;
    nanInds = find( imag( c_gr ) ~= 0 );
    c_gr( nanInds ) = NaN;
    
    % Plot dispersion curves
    cphPlot = plot( fVector, c_ph, 'k', 'LineWidth', 2.2 );
    cgrPlot = plot( fVector, c_gr, '--k', 'LineWidth', 2.2 );
    
end

% Format plot
legHandle = legend( [cphPlot, cgrPlot], ...
    '\,\,Phase Speed $c_{\rm ph}$', ...
    '\,\,Group Speed $c_{\rm gr}$' );
set( legHandle, ...
    'Location', 'SouthEast', ...
    'Interpreter', 'latex', ...
    'FontSize', 24);
xlabel( 'Frequency [Hz]', 'FontSize', 20 );
ylabel( 'Sound Speed [m/s]', 'FontSize', 20 );
ylim([0, 2500]);
zoom xon;

% Add text boxes to label modes
annotation( 'textbox', ...
    'string', 'Mode 1', ...
    'interpreter', 'latex', ...
    'EdgeColor', 'none', ...
    'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 16, ...
    'BackgroundColor', 'w' );
annotation( 'textbox', ...
    'string', 'Mode 2', ...
    'interpreter', 'latex', ...
    'EdgeColor', 'none', ...
    'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 16, ...
    'BackgroundColor', 'w' );
annotation( 'textbox', ...
    'string', 'Mode 3', ...
    'interpreter', 'latex', ...
    'EdgeColor', 'none', ...
    'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 16, ...
    'BackgroundColor', 'w' );
annotation( 'textbox', ...
    'string', 'Mode 4', ...
    'interpreter', 'latex', ...
    'EdgeColor', 'none', ...
    'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 16, ...
    'BackgroundColor', 'w' );