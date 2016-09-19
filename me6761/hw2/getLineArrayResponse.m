% Problem 2 (function so we can have subroutines)

function [ arrayOutput, thetaVector ] = ...
    getLineArrayResponse( source, sensors )

% Get values from input structs
tVector = source.tVector;
t0 = source.t0;
tMax = source.tMax;
fVector = source.fVector;
f0 = source.f0;
BW = source.BW;
sourcePosition = source.position;

numSensors = sensors.number;
dtheta = sensors.dtheta;
c0 = sensors.soundSpeed;

% Get sensor positions along axis
if isnan( sensors.positions )
    spacing = sensors.spacing;
    sensorPositions = getSensorPositions( numSensors, spacing );
else
    sensorPositions = sensors.positions;
end

% Get the signal received be each sensor
vectorsToSource = bsxfun( @minus, sensorPositions, sourcePosition );
sensorTimeDelays = sqrt( sum( vectorsToSource.^(2) ) )./c0;
receivedData = zeros( numSensors, length(tVector) ); % Initialize
for sensorCount = 1:numSensors
    currentT = tVector - t0 - sensorTimeDelays( sensorCount );
    receivedData( sensorCount, : ) = gauspuls( currentT, f0, BW );
end

%%%%%%%% DEBUG %%%%%%%%%%%%
figure()
hold all;
interval = 2.5;
for count = 1:numSensors
    plot( tVector, receivedData(count, :) + interval.*(count - 1) );
end
xlabel( 'Time [s]' );
set( gca, ...
    'YTick', 0:interval:numSensors.*interval, ...
    'YTickLabel', '' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the frequency domain signal
receivedDataTilde = fft( receivedData );

% Compute angle vector
thetaVector = 0:dtheta:pi;

% Compute the beamformer output for each look angle
for thetaCount = 1:length(thetaVector)
    
    % Get look direction unit vector
    theta = thetaVector(thetaCount);
    lookDirection = [ cos(theta); sin(theta); 0 ];
    
    % Compute the associated delays
    tau = sum(bsxfun(@times, lookDirection', sensorPositions'), 3)./c0;
    tau = tau(:, 1); % Get rid of null columns
    
    % Now that we have the time delays, apply them to the source by
    % way of a factor of e^(j*omega*tau) at each frequency
    for freqCount = 1 : floor(length(fVector)./2)
        factor = exp( 2.*pi.*fVector(freqCount).*tau );
        % Apply the factor to the output
        arrayResponse( freqCount, thetaCount ) = ...
            abs( factor'*receivedDataTilde( :, freqCount ) );
    end
end

% Convert back to time domain
arrayOutput = ( 2./length(tVector) )*real( ifft(arrayResponse, [], 1) );

end

% --------------------- SUBROUTINES --------------------------
function [sensorPositions] = getSensorPositions( numSensors, spacing )

% If the number of positions is even, place the origin between the two
% middle sensors. Otherwise, place the origin at the middle sensor.
numSensorsEven = ( numSensors./2 - floor(numSensors./2) == 0 );
if numSensorsEven
    xStart = (-(numSensors-1)./2).*spacing;
    xEnd = ((numSensors-1)./2).*spacing;
    xPos = xStart:spacing:xEnd;
else
    xStart = (-numSensors./2).*spacing;
    xEnd = (numSensors./2).*spacing;
    xPos = xStart:spacing:xEnd;
end

% Sensors are along the x-axis
yPos = 0.*xPos;
zPos = 0.*xPos;
sensorPositions = [ xPos; yPos; zPos ];

end