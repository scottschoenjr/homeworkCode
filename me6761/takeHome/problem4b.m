%**************************************************************************
%
% ME-6761 Final Project, Problem 4b
%
%   Propagation of pulse in 2D waveguide -- Individual Modes
%
%     Scott Schoen Jr 20161130
%
%**************************************************************************

clear all
close all
clc

% Channel parameters
H = 30; % Depth [m]
c = 1500; % Sound speed [m/s] 
tStart = 0; % [s]
tEnd = 40; % [s]

% Source and receiver
xSrc = 0; % [m]
ySrc = 8; % [m]
xRec = 25E3; % [m]
yRec = 21; % [m]

% Pulse paramters
fc = 200; % Center frequency [Hz]
tau = 0.0125; % Width parameter [s]
t0 = 35E-3; % Delay [s]

% Reconstruction parameters (Used to save time)
Fs = 33.*pi.*fc; % [Hz]
fMin = 0; % [Hz]
fMax = Fs; % [Hz]
nModes = 4;
evanescentModesToKeep = 0;

% Compute parameters
dt = 1./Fs; % Time step [s]
tVector = tStart : dt : tEnd; % Time vector [s]
fVector = linspace( 0, Fs, length(tVector) ); % Frequency vector [Hz]
omegac = 2.*pi.*fc; % Angular frequency [rad/s]
rCurs = sqrt( (xSrc - xRec).^(2) + (ySrc - yRec).^(2) ); % Distance [m]
r = abs(xSrc - xRec); % Horizontal range to receiver [m]

% Create pulse
s = sin( omegac.*(tVector - t0) ).*exp( -( (tVector - t0)./tau ).^(2) );

% Take FFT of pulse
sTilde = fft(s);

% Initialize received signal/FFT of received signal
sRecTilde = 0.*sTilde;

% Initialize individual mode matrix
modalContributions = zeros( nModes, length(tVector) );

% Get the contribution at each frequency
tic
for fCount = 1:length(fVector)
    
   % Compute calculations only within bandwidth to save time
   % NOTE - This assumes fVector increases monotonically!
   f = fVector( fCount );
   if f < fMin
       continue;
   elseif f > fMax
       continue; % Break not allowed in parfor
   end
   
   % Initialize Green's function at this frequency
   Gtilde = 0;
   
   % We'll include propagating modes as well as the first few evanescent
   % modes to capture some near field if desired.
   omega = 2.*pi.*f;
   k = omega./c;
   evanescentModesKept = 0;
   
   for modeCount = 1:nModes
       
       % Get propagating mode number
       n = modeCount;
       ky = (n - 0.5).*pi./H;
       kx = sqrt( k.^(2) - ky.^(2) );
       
       % If kx == 0, add tiny amount so we don't blow up
       if kx == 0
           continue;
       end

       % Keep evanescent modes only if we're below the number we want to
       % keep.
       if ~isreal( kx )
           evanescentModesKept = evanescentModesKept + 1;
           if evanescentModesKept > evanescentModesToKeep
               break;
           end
       end
       
       % Add in contribution
       psi  = sin( ky.*yRec );
       psi0 = sin( ky.*ySrc );
       denominator = kx;
       modeContribution = (4./(1j.*H)).*( ...
           psi.*psi0./denominator ...
           ).*exp( 1j.*kx.*r );
       Gtilde = Gtilde + modeContribution;
       
       % Save contribution
       modalContributionsTilde( modeCount, fCount ) = ...
           sTilde(fCount).*modeContribution;
       
   end
   
   % Add this frequency's contribution to the total field
   sRecTilde(fCount) = sTilde(fCount).*Gtilde;
    
end
toc

% Get the received signal in the time domain
sRec = fliplr( ifft( sRecTilde ) );

% Normalize
sNorm = s./max(abs(s));
sRecNorm = sRec./max(abs(real(sRec)));

% Normalize
sNorm = s./max(abs(s));
sRecNorm = sRec./max(max(abs(real(sRec))));

% Plot of modal contributions at receiver
figure()
grayColor = 0.6.*[1, 1, 1];
hold all;
plotOffset = 2.5;
offsets = 0 : plotOffset : nModes.*plotOffset;
for modeCount = 1:nModes
    % Stack on top of each other
    offset = offsets( modeCount );
    recData = modalContributionsTilde( modeCount, : );
    yData = fliplr( ifft( recData ) );
    yDataNorm = real(yData)./max(abs(yData)) + offset;
    plot( tVector, yDataNorm, 'Color', grayColor);
end
% Plot total received signal for reference
offset = offsets( end );
plot( tVector, 1.5.*real(sRecNorm) + offset, 'k' );

set( gca, 'YTick', offsets, ...
    'YTickLabels', {'1';'2';'3';'4'; 'Total'});
ylabel( 'Mode Number', 'FontSize', 20 );

xlabel( 'Time [s]', 'FontSize', 22);
xlim([10, 20]);

title( sprintf( 'Range = %d km', round(xRec./1E3) ), ...
    'FontSize', 22);

box on;
zoom xon;

xlim( [10, 20] );