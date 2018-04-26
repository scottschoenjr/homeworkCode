% *************************************************************************
% Source Pulse
%
%   Function to define Gaussian-windowed pulse (similar to gauspuls) and 
%   its second derivative.
%
% Inputs
%   tVec - Time vector [s]
%   P0   - Pressure amplitude [Pa]
%   f0   - Center frequency [Hz]
%   t0   - Time shift [s] [optional]
%   tw   - Time width [s] [optional]
% 
% Ouputs
%   s    - Signal [Pa]
%   s_tt - Second derivative of signal [Pa/s^2]
%
% *************************************************************************

function [ s, s_tt ] = pulseWaveform( tVec, P0, f0, t0, tw )

% Set defaults
if nargin < 6
    tw = 3/f0; % Pulse width [s]
end
if nargin < 5
    t0 = 6/f0; % Pulse offset [s]
end

omega0 = 2.*pi.*f0; % Angular frequency

% Created shifted time
t = tVec - t0;

% Define time series
s = P0.*exp(-(t./(tw/2)).^2).* sin(omega0.*t);

% Define its second derivative
s_tt = ...
    -8.*P0.*exp( -4.*t.^(2)./tw.^(2) ).*sin(omega0.*t)./tw.^(2) + ...
    64.*P0.*t.^(2).* exp( -4.*t.^2./tw.^2 ).*sin( omega0.*t )./tw.^(4) - ...
    32.*P0.*t.*exp(-4. *t.^2./tw.^2).*cos( omega0.*t ).*pi.*f0./tw.^(2) - ...
    4.*P0 .* exp(-4.*(tVec-t0).^(2)./tw.^2) .* sin(omega0.*t).*pi.^(2).*f0.^(2);

% Force excitation to 0 at end
zeroIndex = find(tVec > 12./f0, 1);
s( zeroIndex : end ) = 0;
s_tt( zeroIndex : end ) = 0;

end

