% Example from class

clear all
close all
clc

% Velocity

% Declare symbolic variables
syms vA betadot psidot 

% Define the vector connecting the ball and socket to the collar
rB = [ 1; 2; 0 ]; % Collar [m]
rA = [ 0; 0; 2 ]; % Ball and socket
r_AB = rA - rB;

% Define the unit vectors
e_x = [1; 0; 0];
e_y = [0; 1; 0];
e_z = [0; 0; 1];

% Define the velocity along the y axis of the collar
vB = 3*e_y;

% Define the unit vector of the pin around which the fork rotates. Since it
% must be perpendicular to both AB and y, we can just cross these two. To
% have it point in the positive x(ish) direction, put y first.
e_pin = cross( e_y, r_AB );
e_pin = e_pin./( norm( e_pin ) ); % Normalize it

% Now define the total rotation. We have beta about the pin axis, but also
% psi about the y axis. Thus
omega = betadot*e_pin + psidot*e_y;

% Now since the velocity of the two ends must match, get all terms on one
% side and use the equation solver
velocityEquation = vB + cross( omega, r_AB ) + vA*e_z;

% Use the solver to solve for the desired variables
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), velocityEquation(3), ...
    vA, betadot, psidot ... % Variables to solve for
    );

vA_computed = double( velocitySolutionStruct.vA );
betadot_computed = double( velocitySolutionStruct.betadot );
psidot_computed = double( velocitySolutionStruct.psidot );

% Print results
fprintf( 'v_A      = %0.4f m/s\n', vA_computed );
fprintf( 'thetadot = %0.4f rad/s\n', betadot_computed );
fprintf( 'psidot   = %0.4f rad/s\n\n', psidot_computed );

% Acceleration

% Declare symbolic unknowns
syms vadot betaddot psiddot

% Define the solved-for angluar velocity
omegaAB = betadot_computed*e_pin + psidot_computed*e_y;

% Define the angular acceleration alpha = d/dt( omega )
alpha = betaddot*e_pin + betadot_computed*cross( omegaAB, e_pin ) + psiddot*e_y;

% Since the velocity of the collar B is constant set aB to 0
aB = 0*e_z;

% Now set up the full acceleration expression, move the B terms to the
% right, and use the MATLAB solver
accelerationEquation = aB + cross( alpha, r_AB ) ...
    + cross( omegaAB, cross( omegaAB, r_AB ) ) + vadot*e_z;

accelerationSolutionStruct = solve( ...
    accelerationEquation(1), accelerationEquation(2), accelerationEquation(3), ...
    vadot, betaddot, psiddot ... % Variables to solve for
    );

vadot = double( accelerationSolutionStruct.vadot );
betaddot = double( accelerationSolutionStruct.betaddot );
psiddot = double( accelerationSolutionStruct.psiddot );

% Print results
fprintf( 'a_A       = %0.4f m/s^2\n', vadot );
fprintf( 'thetaddot = %0.4f rad/s^2\n', betaddot );
fprintf( 'psiddot   = %0.4f rad/s^2\n', psiddot );






