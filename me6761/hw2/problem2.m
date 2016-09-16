% Problem 2 (function so we can have subroutines)

function [] = problem2()

clear all
close all
clc

% Define array positions along the x-axis
numSensors = 16; % Number of sensors
spacing = 15;    % Spacing between sensors [m]
dtheta = 1.*pi./180; % Angular resolution [rad]
c0 = 1500; % Sound Speed [m/s]

% Get sensor positions along axis
sensorPositions = getSensorPositions( numSensors, spacing );

% Compute angle vector
thetaVector = 0:dtheta:pi;

% Compute the beamformer output for each look angle
for count = 1:length(thetaVector)
    
    % Get look direction unit vector
    theta = thetaVector(count);
    lookDirection = [ cos(theta); sin(theta); 0 ]; 
    
    % Dot each sensor position with the look direction to get the
    % associated time delays
    directionalDelays = ...
        bsxfun(@times, lookDirection, sensorPositions)./c0;
    % The total delay is the sum of the delays in each direction
    timeDelays = sum( directionalDelays', 2);

end

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

% Sensors are in the x-y plane
yPos = 0.*xPos;
zPos = 0.*xPos;
sensorPositions = [ xPos; yPos; zPos ];

end