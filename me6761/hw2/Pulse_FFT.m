% Problem 1(c)

clear all
close all
clc

% Set parameters
Fs = 500; % [Hz]
fc = 50; % [Hz]
T = 100; % [s]
t0 = 50; % [s]

% Set tau values
tauValues = [ 10, 1, 0.1 ]; % [s]

% Define function for each tau
omegac = 2.*pi.*fc;
dt = 1./Fs;
tVector = 0:dt:T;
fVector = linspace( 0, Fs, length(tVector) );

% Compute FFT for each
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
        '~$\tau$ = ', num2str(tau) ];
    
end

% Format plot
subplot(2,1,1)
ylabel('Normalized Magnitude');
set(gca, 'XTickLabel', '' );
set(gca, 'Position', [0.1, 0.47, 0.8, 0.3349] );
xlim( [0, 250] );
legend(legString, 'interpreter', 'latex');
box on;

subplot(2, 1, 2)
ylabel( 'Phase [°]' );
set( gca, 'YTick', -180:90:180 );
xlabel('Frequency [Hz]');
xlim( [0, 250] );
box on;
