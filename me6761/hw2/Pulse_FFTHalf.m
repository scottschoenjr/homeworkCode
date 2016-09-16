% Problem 1 (part e)

clear all
close all
clc

% Set parameters
Fs = 500; % [Hz]
fc = 50; % [Hz]
T = 100; % [s]
t0 = 50; % [s]
omegac = 2.*pi.*fc; % [rad/s]
dt = 1./Fs; % [s]
tVector = 0:dt:T; % [s]
tau = 0.1; % [s]

% Define signal
arg = ( (tVector - t0)./tau ).^(2);
signal = sin( omegac.*tVector ).*exp( -arg );

% Compute FFT
X = fft( signal );

% Compute first half of signal
halfIndex = floor( length(X)./2 );
XHalf = 0.*X;
XHalf(1:halfIndex) = X(1:halfIndex);

% Compute reconstructed signal
Fhalf = 2.*real( ifft( XHalf ) );

% Plot results
figure(15);
hold all;

% Plot signal
plot( tVector, signal - Fhalf );
% plot( tVector, Fhalf );

% Format plot
xlabel('Time [s]');
xlim( [40, 60] );
ylabel('$f - f_{\rm half}$ [AU]', 'interpreter', 'latex', 'FontSize', 16);
% legend( ' Original', ' Reconstructed');
box on;
