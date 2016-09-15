% Problem 1 (part f)

clear all
close all
clc

% Set parameters
Fs = 500; % [Hz]
fc = 50; % [Hz]
T = 100; % [s]
t0 = 50; % [s]
tau = 0.1;

% Create figure
figure(161);
hold all;

% Set tau values
tauValues = [ 10, 1, 0.1 ]; % [s]

% Define function for each tau
omegac = 2.*pi.*fc;
dt = 1./Fs;

% Define signal
arg = ( (tVector - t0)./tau ).^(2);
signal = sin( omegac.*tVector ).*exp( -arg );

% Plot signal
plot( tVector, signal );

% Define the frequency vector
fVector = linspace( 0, Fs, length(tVector) );

% Compute and plot the FFT for each value of tau
t0 = 50; % Reset t0 [s]

figure(162)
set(gcf, 'Position', [100, 50, 1000, 850] );
    
    % Take FFT of signal
    signalTilde = fft( signal );
    
    % Plot signal magnitude and phase
    subplot(2, 1, 1);
    hold all;
    plot( fVector, abs(signalTilde)./max(abs(signalTilde)) );
    ylabel('Normalized Magnitude');
    set(gca, 'XTickLabel', '' );
    xlim( [0, 250] );
    
    subplot(2, 1, 2);
    hold all;
    plot( fVector, 180.*angle(signalTilde)./pi );
    
    % Add legend entry
    legString{count} = [ ...
        '~$\tau$ = ', num2str(tau) ];
    

legend(legString, 'interpreter', 'latex');
box on;

subplot(2, 1, 2)
ylabel( 'Phase [�]' );
set( gca, 'YTick', -180:90:180 );
xlabel('Frequency [Hz]');
xlim( [0, 250] );
box on;
