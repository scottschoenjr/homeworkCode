% Problem 3d - Change Frequency of second source

clear all;
close all;
clc;

% -----------------------------------------------------------------------

% Define time and signal parameters
Fs = 400;   % Sampling frequency [Hz]
tMax = 10;   % Time to record for [s]
t0 = 1;      % Time of the center of the signal [s]
f01 = 50;    % Center frequency of the first signal [Hz]
f02 = 25;    % Center frequency of the first signal [Hz]
BW = 0.1;    % Fractional bandwidth of the signal
range = 7E3; % Range of the signal relative to the array center [m]
bearing1 = 45.*pi./180; % Bearing of the source re array center [rad]
bearing2 = 38.7.*pi./180; % Bearing of the source re array center [rad]

% Define sensor parameters
numSensors = 16;     % Number of sensors
spacing = 30;        % Spacing between sensors [m]
dtheta = 0.5.*pi./180; % Angular resolution [rad]
c0 = 1500;           % Sound Speed [m/s]

% -------------------------------------------------------------------------

% Assemble first source struct
dt = 1./Fs;
source1.tVector = 0:dt:tMax;
source1.t0 = t0;
source1.tMax = tMax;
source1.BW = BW; % 0 will give plane wave
source1.fVector = linspace( 0, Fs, length(source1.tVector) );
source1.f0 = f01;
source1.amplitude = 1;
source1.position = [ range.*sin( bearing1 ); range.*cos( bearing1 ); 0 ];
% And second source struct
dt = 1./Fs;
source2.tVector = 0:dt:tMax;
source2.t0 = t0;
source2.tMax = tMax;
source2.BW = BW; % 0 will give plane wave
source2.fVector = linspace( 0, Fs, length(source2.tVector) );
source2.f0 = f02;
source2.amplitude = 1;
source2.position = [ range.*sin( bearing2 ); range.*cos( bearing2 ); 0 ];
% Specify Noise Source struct
dt = 1./Fs;
noise.tVector = 0:dt:tMax;
noise.t0 = t0;
noise.tMax = tMax;
noise.BW = -1; % Negative value to provide noise
noise.fVector = linspace( 0, Fs, length(source1.tVector) );
noise.f0 = f01;
noise.position = [ range.*sin( bearing1 ); range.*cos( bearing1 ); 0 ];
noise.amplitude = 3.16E-4; % = 50dB [Pa]

% Assemble sensor struct
sensors.number = numSensors;
sensors.spacing = spacing;
sensors.positions = NaN; % Specify non-uniform or other sensor spacing here
sensors.dtheta = dtheta;
sensors.soundSpeed = c0;

% Get the array response as a function of theta
allSources = [source1, source2, noise];
[ arrayResponse1, thetaVector1 ] = ...
    getLineArrayResponse( [source1, noise], sensors );
[ arrayResponse2, thetaVector2 ] = ...
    getLineArrayResponse( [source2, noise], sensors );
[ arrayResponseTotal, thetaVectorTotal ] = ...
    getLineArrayResponse( allSources, sensors );

% Plot results
figure()
hold all
totalPlot = plot( 180.*thetaVectorTotal./pi, ...
    max(arrayResponseTotal)./max(max(arrayResponseTotal)), 'k' );
source1Plot = plot( 180.*thetaVector1./pi, ...
    max(arrayResponse1)./max(max(arrayResponseTotal)), '--k', 'LineWidth', 1 );
source2Plot = plot( 180.*thetaVector2./pi, ...
    max(arrayResponse2)./max(max(arrayResponseTotal)), ':k', 'LineWidth', 1 );
xlabel( 'Look Angle \theta [deg]' );
xlim( [-90, 90] );
set(gca, 'XTick', -90:30:90 );
ylabel( 'Normalized Response' );
ylim([0, 1.2]);
% title( sprintf( 'f_{0} = %0.2f Hz', f01 ) );
legend( [source1Plot, source2Plot, totalPlot], ...
    ' Source 1', ' Source 2', ' Total', ...
    'Location', 'NorthWest');
box on;