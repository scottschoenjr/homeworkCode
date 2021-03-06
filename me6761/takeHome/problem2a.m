%**************************************************************************
%
% ME-6761 Final Project, Problem 2
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

% Set range of tau values to sweep over
numTaus = 3;
tauVector = linspace( 0.030, 0.040, numTaus );

% Create time vector for each pulse
tVector = linspace( 0, 1, 50E3);
dt = tVector(2) - tVector(1);
fVector = linspace( 0, 1./dt, length(tVector) );

% Set lineshapes
lineShapes = { '--k'; '-k'; '-.k' }; 

% Compute the pulse for each tau
for tauCount = 1:numTaus
    
    tau = tauVector(tauCount);
    f = sin( 2.*pi.*fc.*tVector ).*exp( -(tVector./tau).^(2) );
    F = fft(f);
    
    % Plot spectrum
    figure(1)
    hold all;
    plot( fVector, 20.*log10(abs(F)./max(abs(F))), lineShapes{tauCount} )
    legendString{ tauCount } = ...
        [' $\tau =$ ', sprintf('%0.1f', 1E3.*tau), ' ms' ];
    
end

% Format plot
box on;
ylim( [-6, 0] );
xlim( [fc - 18, fc + 18] );
zoom xon

xlabel('Frequency [Hz]');
ylabel('Normalized Amplitude [dB]');
legend( legendString, 'interpreter', 'latex' );