%**************************************************************************
%
% ME-6761 Final Project, Problem 2b
%
%   Plotting Gaussian-windowed pulse
%
%     Scott Schoen Jr 20161114
%
%**************************************************************************

clear all
close all
clc

% Pulse paramters
fc = 100; % Center frequency [Hz]
tau = 0.0125; % [s]

% Set range of tau values to sweep over
numt0s = 500;
t0Vector = linspace( 0, 50E-3, numt0s );

% Create time vector
tVector = linspace( 0, 1, 10E3);

% Initialize vector to hold value at 0
valueAt0 = zeros(1, numt0s);

% For each delay t0
for t0Count = 1:numt0s
    
    % Compute the pulse
    t0 = t0Vector( t0Count );
    f = sin( 2.*pi.*fc.*(tVector - t0) ).*exp( -( (tVector - t0)./tau ).^(2) );
    
    % Normalize and get value at 0
    fNorm = f./max(abs(f));
    valueAt0(t0Count) = fNorm(1);
    
%     %%%%%%% DEBUG %%%%%%%
%     figure(999)
%     hold all;
%     plot( tVector, fNorm );
%     xlabel('Time [s]');
%     ylabel('Amplitude' );
%     %%%%%%%%%%%%%%%%%%%%%
      
end

% Format plot
figure()
hold all;
box on;

% Set shaded region
xMin = 30;  % [ms]
xSpan = 10; % [ms]

rectangle( ...
    'Position', [ xMin, -120, xSpan, 240 ], ...
    'FaceColor', 0.8.*[1, 1, 1], ...
    'EdgeColor', 'none' ...
    );

plot( 1E3.*t0Vector, 100.*valueAt0, 'k' );
xlabel('$t_{0}$ [ms]');
ylabel('Relative Value at $t = 0$ [\%]');

ylim( [-100, 100] );
xlim( [0, 50] );

% Plot inset
figure()
box on;
hold all;

plot( 1E3.*t0Vector, 100.*valueAt0, 'k' );

% Plot limit lines
plot( [xMin, xMin + xSpan], [ 0.1, 0.1 ], '--k', 'LineWidth', 0.5 );
plot( [xMin, xMin + xSpan], [ -0.1, -0.1 ], '--k', 'LineWidth', 1 );

xlabel('$t_{0}$ [ms]');
ylabel('Relative Value at $t = 0$ [\%]');

ylim( [-0.2, 0.2] );
xlim( [xMin, xMin + xSpan] );