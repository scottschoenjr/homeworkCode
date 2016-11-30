%**************************************************************************
%
% ME-6761 Final Project, Problem 4
%
%   Propagaition of pulse in 2D waveguide
%
%     Scott Schoen Jr 20161129
%
%**************************************************************************

clear all
close all
clc

% Channel parameters
H = 30;     % Depth [m]
c = 1500;   % Sound speed [m/s]
tStart = 0; % [s]
tEnd = 50;  % [s]

% Source position
xSrc = 0; % [m]
ySrc = 8; % [m]

% Receiver positions
dx = 100;
xRecVector = 20E3 : dx : 30E3;
numReceivers = length( xRecVector );
yRecVector = 21.*ones( 1, numReceivers );

% Pulse paramters
fc = 200; % Center frequency [Hz]
tau = 0.0125; % Width parameter [s]
t0 = 35E-3; % Delay [s]

% Reconstruction parameters (Used to save time)
fMin = 0; % [Hz]
fMax = 20E3; % [Hz]
nModes = 4;
evanescentModesToKeep = 0;

% Compute parameters
Fs = 3.*pi.*fc; % Sampling rate [Hz]
dt = 1./Fs;      % Time step [s]
tVector = tStart : dt : tEnd;                 % Time vector [s]
fVector = linspace( 0, Fs, length(tVector) ); % Frequency vector [Hz]
omegac = 2.*pi.*fc; % Angular frequency [rad/s]

% Create pulse
s = sin( omegac.*(tVector - t0) ).*exp( -( (tVector - t0)./tau ).^(2) );

% Take FFT of pulse
sTilde = fft(s);

% Initialize received signal/FFT of received signal
sRecTilde = 0.*sTilde;
sRec = zeros( numReceivers, length(tVector) );

tic;
% Compute the received signal at each receiver. This could be slightly
% quicker, if we moved this inside the loop over mode numbers, but speed
% is less of a concern at the moment.
for recCount = 1:numReceivers
        
    % Get distance values for this receiver
    xRec = xRecVector( recCount );
    yRec = yRecVector( recCount );
    rCurs = sqrt( (xSrc - xRec).^(2) + (ySrc - yRec).^(2) );
    r = abs(xSrc - xRec); % Horizontal range to receiver [m]
      
    % Get the contribution at each frequency
    for fCount = 1:length(fVector)
        
        % Compute calculations only within bandwidth to save time
        % NOTE - This assumes fVector increases monotonically!
        f = fVector( fCount );
        if f < fMin
            continue;
        elseif f > fMax
           continue;
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
            
            % Keep evanescent modes only if we're below the number we
            % want to keep.
            if ~isreal( kx )
                evanescentModesKept = evanescentModesKept + 1;
                if evanescentModesKept > evanescentModesToKeep
                    break;
                end
            end
            
            % Compute mode shapes for this frequency and position
            psi  = sin( ky.*yRec );
            psi0 = sin( ky.*ySrc );
            denominator = kx;
            modeContribution = (4./(1j.*H)).*( ...
                psi.*psi0./denominator ...
                ).*exp( 1j.*kx.*r );
            Gtilde = Gtilde + modeContribution;
            
        end
        
        % Add this frequency's contribution to the total field
        sRecTilde(fCount) = sTilde(fCount).*Gtilde;
        
    end
    
    % Get the received signal in the time domain
    sRec( recCount, : ) = fliplr( ifft( sRecTilde ) );
    
end
toc

% Normalize
sNorm = s./max(abs(s));
sRecNorm = sRec./max(max(abs(real(sRec))));

%% Format plots

% Plot of all receivers' signals
figure()
hold all;
box on;

% Plot signal
[xPlot, tPlot] = meshgrid( xRecVector, tVector );

% Get large variables out of memory
clear s sTilde sNorm sRecTilde fVector

% Plot received signal
rangePlot = pcolor( tPlot', xPlot'./1E3, abs(hilbert(real(sRecNorm))) );
xlabel( 'Time [s]' );
ylabel('Recever Position [km]');

xlim( [12, 22] );
ylim( [20, 30] );

blackAndWhite = flipud( colormap( gray ) );
colormap( blackAndWhite );
set(rangePlot, 'EdgeColor', 'none' );

box on;
zoom xon;




