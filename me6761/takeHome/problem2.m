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
numTaus = 5;
tauVector = linspace( 0.010, 0.015, numTaus );

% Create time vector for each pulse
tVector = linspace( 0, 1, 20000);
dt = tVector(2) - tVector(1);
fVector = linspace( 0, 1./dt, length(tVector) );

% Compute the pulse for each tau
for tauCount = 1:numTaus
    
    tau = tauVector(tauCount);
    f = sin( 2.*pi.*fc.*tVector ).*exp( -tVector./tau );
    F = fft(f);
    
    % Plot spectrum
    figure(1)
    hold all;
    plot( fVector, 20.*log10(abs(F)./max(abs(F))) )
    legendString{ tauCount } = [' \tau = ', sprintf('%0.3f', tau) ];
    
end

% Format plot
box on;
ylim( [-3, 0] );
xlim( [fc - 15, fc + 15] );
zoom xon

xlabel('Frequency [Hz]');
ylabel('Normalized Amplitude [dB]');
legend( legendString );