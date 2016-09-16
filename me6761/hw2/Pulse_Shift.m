% Problem 1 (part f)

clear all
close all
clc

% Set parameters
Fs = 500; % [Hz]
fc = 50; % [Hz]
T = 100; % [s]
t0 = 50; % [s]
tau = 0.1; % [s]
tShift = 10; % [s]
omegac = 2.*pi.*fc; % [rad/s]
dt = 1./Fs; % [s]
tVector = 0:dt:T; % [s]
fVector = linspace( 0, Fs, length(tVector) ); % [Hz]
omegaVector = 2.*pi.*fVector; % [rad/s]

% Define signal and shifted signal
arg = ( (tVector - t0)./tau ).^(2);
signal = sin( omegac.*tVector ).*exp( -arg );
argShifted = ( (tVector - t0 - tShift)./tau ).^(2);
signalShifted = sin( omegac.*(tVector - tShift) ).*exp( -argShifted );

% Compute the fft of the signal
X = fft( signal );
Xshifted = X.*exp(-1j.*omegaVector.*tShift);

% Reconstruct the signal with theime delay
signalRecon = ifft( X );
signalShiftedRecon = ifft( Xshifted );

% Create figure
figure(161);
hold all;

% Plot signal
plot( tVector, signalShifted );
plot( tVector, real(signalShiftedRecon) );

xlabel('Time [s]');
xlim( [50, 70] );
ylabel('$x(t-\tau) - F^{-1}[e^{-j\omega\tau}X(\omega)]$ [AU]', 'interpreter', 'latex', 'FontSize', 16);
% legend( ' Original', ' Reconstructed');
box on;


