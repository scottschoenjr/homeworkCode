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
pulseSignal = ( BW ~= 0 );

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
    
    if pulseSignal
        receivedData( sensorCount, : ) = gauspuls( currentT, f0, BW );
    else
        receivedData( sensorCount, : ) = cos( 2.*pi.*f0.*currentT );
    end
end

% Get the frequency domain signal
receivedDataTilde = fft( receivedData' );
receivedDataTilde = receivedDataTilde'; % Keep dimensions of receivedData

% Compute angle vector note that this is polar angle
thetaVector = -pi:dtheta:pi;

% Compute the beamformer output for each look angle
for thetaCount = 1:length(thetaVector)
    
    % Get look direction unit vector
    theta = thetaVector(thetaCount);
    lookDirection = [ sin(theta); cos(theta); 0 ];
    
    % Compute the associated delays
    tau = sum(bsxfun(@times, lookDirection', sensorPositions'), 3)./c0;
    tau = tau(:, 1); % Get rid of null columns
    
    % Now that we have the time delays, apply them to the source by
    % way of a factor of e^(j*omega*tau) at each frequency
    for freqCount = 1 : floor(length(fVector)./2)
        factor = exp( 1j.*2.*pi.*fVector(freqCount).*tau );
        % Apply the factor to the output
        arrayResponse( freqCount, thetaCount ) = ...
            abs( factor'*receivedDataTilde( :, freqCount ) );
    end
end

% Convert back to time domain (why fiip?)
arrayOutput = fliplr(( 2./length(tVector) )*real( ifft(arrayResponse, [], 1) ));

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