%**************************************************************************
%
% ME-6761 Final Project, Problem 3
%
%   Propagaition of pulse in 2D waveguide
%
%     Scott Schoen Jr 20161116
%
%**************************************************************************

clear all
close all
clc

% Channel parameters
H = 30; % Depth [m]
c = 1500; % Sound speed [m/s] 
tStart = 0; % [s]
tEnd = 20; % [s]
Fs = 20E3; % [Hz]

% Source and receiver
xSrc = 0; % [m]
ySrc = 8; % [m]
xRec = 20E3; % [m]
yRec = 29; % [m]

% Pulse paramters
fc = 100; % Center frequency [Hz]
tau = 0.0125; % Width parameter [s]
t0 = 35E-3; % Delay [s]

%%%%%%% DEBUG %%%%%%%
fc = 1E3; % Center frequency [Hz]
tau = 2E-3; % Width parameter [s]
t0 = 35E-3; % Delay [s]
%%%%%%%%%%%%%%%%%%%%%%

% Reconstruction parameters (Used to save time)
fMin = 0; % [Hz]
fMax = 20E3; % [Hz]
nModes = 50;
evanescentModesToKeep = 20;

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
       denominator = k.^(2) - kx.^(2);
       Gtilde = Gtilde - (4.*pi./H).*( ...
           psi.*psi0./denominator ...
           ).*exp( 1j.*kx.*r );
       
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

% Format plot
figure()
hold all;
box on;

plot( tVector, sNorm + 2, '--k' );
plot( tVector, real(sRecNorm), 'k' );
xlabel('$t$ [s]');
ylabel('$s(t)$');

ylim( [-2, 5] );
xlim( [0, 20] );

