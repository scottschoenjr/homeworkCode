%**************************************************************************
%
% ME-6761 Final Project, Problem 4c
%
%   Computing phase and group speeds for first four modes
%
%     Scott Schoen Jr 20161130
%
%**************************************************************************

clear all
close all
clc

% Channel paramters
H = 30; % [m]
f = 200; % [Hz]
c = 1500; % [m/s]

% Compute first few mode shapes
numModes = 4;
omega = 2.*pi.*f;
k = omega./c;

% Initialize
c_gr = zeros( 1, numModes );
c_ph = zeros( 1, numModes );

% Plot each mode
for modeCount = 1:numModes

    % Compute wavenumber
    ky = (modeCount - 0.5).*pi./H;
    
    % Compute phase speed
    c_gr(modeCount) = c.*(1 - (ky./k).^(2)).^(0.5);
    
    % Compute phase speed
    c_ph(modeCount) = c.^(2)./c_gr(modeCount);
    
end

% Print
c_gr
c_ph