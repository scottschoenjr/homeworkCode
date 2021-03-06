% Problem 1(a)

clear all
close all
clc

% Set parameters
Fs = 500; % [Hz]
fc = 50; % [Hz]
T = 100; % [s]
t0 = 50; % [s]

% Part (a) ---------------------

% Create figure
figure(111);
hold all;

% Set tau values
tauValues = [ 10, 1, 0.1 ]; % [s]

% Define function for each tau
omegac = 2.*pi.*fc;
dt = 1./Fs;
tVector = 0:dt:T;
for count = 1:length( tauValues )
   
    tau = tauValues(count);
    arg = ( (tVector - t0)./tau ).^(2);
    
    % Define signal
    signal = sin( omegac.*tVector ).*exp( -arg );
    
    % Plot signal
    plot( tVector, signal );
    
    % Add legend entry
    legString{count} = [ ...
        '~$\tau$ = ', num2str(tau), '~s' ];
    
end

% Format plot
title( ['$t_{0} =~$', num2str(t0), '~s' ], 'interpreter', 'latex' );
xlabel('Time [s]');
xlim( [40, 60] );
ylabel('Signal [AU]');
legend(legString, 'interpreter', 'latex');
box on;

% Shift to t0 = 30 s
t0 = 30; % [s]
figure(112)
hold all;
for count = 1:length( tauValues )
   
    tau = tauValues(count);
    arg = ( (tVector - t0)./tau ).^(2);
    
    % Define signal
    signal = sin( omegac.*tVector ).*exp( -arg );
    
    % Plot signal
    plot( tVector, signal );
    
    % Add legend entry
    legString{count} = [ ...
        '~$\tau$ = ', num2str(tau), '~s' ];
    
end

% Format plot
title( ['$t_{0} =~$', num2str(t0), '~s' ], 'interpreter', 'latex' );
xlabel('Time [s]');
xlim( [20, 40] );
ylabel('Signal [AU]');
legend(legString, 'interpreter', 'latex');
box on;


% Part (c) -----------------------

% Define the frequency vector
fVector = linspace( 0, Fs, length(tVector) );

% Compute and plot the FFT for each value of tau
t0 = 50; % Reset t0 [s]

figure(13)
set(gcf, 'Position', [100, 50, 1000, 850] );
for count = 1:length( tauValues )
   
    tau = tauValues(count);
    arg = ( (tVector - t0)./tau ).^(2);
    
    % Define signal
    signal = sin( omegac.*tVector ).*exp( -arg );
    
    % Take FFT of signal
    signalTilde = fft( signal );
    
    % Plot signal magnitude and phase
    subplot(2, 1, 1);
    hold all;
    plot( fVector, abs(signalTilde)./max(abs(signalTilde)) );
    subplot(2, 1, 2);
    hold all;
    plot( fVector, 180.*angle(signalTilde)./pi );
    
    % Add legend entry
    legString{count} = [ ...
        '~$\tau$ = ', num2str(tau), '~s' ];
    
end

% Format plot
subplot(2,1,1)
ylabel('Normalized Magnitude');
set(gca, 'XTickLabel', '' );
xlim( [0, 250] );
legend(legString, 'interpreter', 'latex');
box on;

subplot(2, 1, 2)
ylabel( 'Phase [�]' );
set( gca, 'YTick', -180:90:180 );
xlabel('Frequency [Hz]');
xlim( [0, 250] );
box on;
