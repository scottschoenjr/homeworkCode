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
tEnd = 15; % [s]
Fs = 20E3; % [Hz]

% Source and receiver
xSrc = 0; % [m]
ySrc = 8; % [m]
xRec = 20E3; % [m]
yRec = 21; % [m]

% Pulse paramters
fc = 100; % Center frequency [Hz]
tau = 0.0125; % Width parameter [s]
t0 = 35E-3; % Delay [s]

% Reconstruction parameters (Used to save time)
fMin = 0; % [Hz]
fMax = 500; % [Hz]
evanescentModesToKeep = 2;

% Compute parameters
dt = 1./Fs; % Time step [s]
tVector = tStart : dt : tEnd; % Time vector [s]
fVector = linspace( 0, Fs, length(tVector) ); % Frequency vector [Hz]
omegac = 2.*pi.*fc; % Angular frequency [rad/s]

% Create pulse
s = sin( omegac.*(tVector - t0) ).*exp( -( (tVector - t0)./tau ).^(2) );

% Take FFT of pulse
sTilde = fft(s);

% Initialize received signal/FFT of received signal
sRecTilde = 0.*sTilde;

% Get the contribution at each frequency
for fCount = 1:length(fVector)
    
   % Compute calculations only within bandwidth to save time
   % NOTE - This assumes fVector increases monotonically!
   f = fVector( fCount );
   if f < fMin
       continue;
   elseif f > fMax
       break;
   end
   
   % Initialize Green's function at this frequency
   Gtilde = 0;
   
   % We'll include propagating modes as well as the first 2 evanescent
   % modes to capture some near field if desired.
   omega = 2.*pi.*f;
   k = omega./c;
   evanescentModesKept = 0;
   
   for modeCount = 1:20
       
       % Get propagating mode number
       n = modeCount;
       ky = (n - 0.5).*pi./H;
       omega_m = c.*ky; % Modal frequency [rad/s]
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
       psi  = sin( ky.*yRec ).*exp(1j.*kx.*xRec);
       psi0 = sin( ky.*ySrc ).*exp(1j.*kx.*xSrc);
       denominator = k.^(2) - omega_m.^(2)./c.^(2);
       Gtilde = Gtilde - 4.*pi.*( ...
           psi.*psi0./denominator ...
           );
       
   end
   
   % Add this frequency's contribution to the total field
   sRecTilde(fCount) = sTilde(fCount).*Gtilde;   
    
end

% Get the received signal in the time domain
sRec = ifft( sRecTilde );

% Normalize
sNorm = s./max(abs(s));
sRecNorm = sRec./max(abs(real(sRec)));

% Format plot
figure()
hold all;
box on;

plot( tVector, sNorm + 2, '--k' );
plot( tVector, real(sRecNorm), 'k' );
xlabel('$t_{0}$ [ms]');
ylabel('$s(t)$');

ylim( [-2, 5] );
xlim( [0, 20] );

