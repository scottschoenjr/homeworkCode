% Problem 2a - Harmonic source

clear all;
close all;
clc;

% -----------------------------------------------------------------------

% Define time and signal parameters
Fs = 400;   % Sampling frequency [Hz]
tMax = 10;   % Time to record for [s]
t0 = 1;      % Time of the center of the signal [s]
f0 = 50;    % Center frequency of the signal [Hz]
BW = 0.1;    % Fractional bandwidth of the signal
range = 7E3; % Range of the signal relative to the array center [m]
bearing = 45.*pi./180; % Bearing of the source re array center [rad]

% Define sensor parameters
numSensors = 16;     % Number of sensors
spacing = 15;        % Spacing between sensors [m]
dtheta = 1.*pi./180; % Angular resolution [rad]
c0 = 1500;           % Sound Speed [m/s]

% -------------------------------------------------------------------------

% Assemble source struct
dt = 1./Fs;
source.tVector = 0:dt:tMax;
source.t0 = t0;
source.tMax = tMax;
source.BW = 0; % 0 will give plane wave
source.fVector = linspace( 0, Fs, length(source.tVector) );
source.f0 = f0;
source.position = [ range.*sin( bearing ); range.*cos( bearing ); 0 ];
source.amplitude = 10; % = 140 dB [Pa]

% Specify Noise Source struct
dt = 1./Fs;
noise.tVector = 0:dt:tMax;
noise.t0 = t0;
noise.tMax = tMax;
noise.BW = -1; % Negative value to provide noise
noise.fVector = linspace( 0, Fs, length(source.tVector) );
noise.f0 = f0;
noise.position = [ range.*sin( bearing ); range.*cos( bearing ); 0 ];
noise.amplitude = 3.16E-4; % = 50dB [Pa]

% Assemble sensor struct
sensors.number = numSensors;
sensors.spacing = 15;
sensors.positions = NaN; % Specify non-uniform or other sensor spacing here
sensors.dtheta = dtheta;
sensors.soundSpeed = c0;

% Get the array response as a function of theta
noiseySource = [source, noise];
[ arrayResponse, thetaVector ] = getLineArrayResponse( source, sensors );

% Plot results
figure()
hold all
plot( 180.*thetaVector./pi, ...
    max(arrayResponse)./max(max(arrayResponse)), 'k' );
xlabel( 'Look Angle \theta [deg]' );
xlim( [-90, 90] );
set(gca, 'XTick', -90:30:90 );
ylabel( 'Normalized Response' );
box on;