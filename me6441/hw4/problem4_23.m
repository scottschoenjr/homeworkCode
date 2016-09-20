% Problem 4-23 in _Engineering Dynamics_ (Ginsberg, 2008)

clear all
close all
clc

% Velocity

% Declare symbolic variables
syms vD betadot psidot

% Define the vector connecting the ball and socket to the collar
rC = [ 0; 3; 0 ];         % Cap C [m]
rD = [ -0.5; 0; 2.5981 ]; % Collar D [m]
r_CD = rC - rD;

% Define the unit vectors
e_X = [1; 0; 0];
e_Y = [0; 1; 0];
e_Z = [0; 0; 1];

% Define rotation speed of the crankshaft
omegaA = (900./60); % [rad/s]
rBC = 1; % [m]

% Define the unit vector for the pin
e_pin = cross( e_Y, r_CD );
e_pin = e_pin./norm(e_pin);

% Define the sides of the velocity equation
LHS = -(rBC*omegaA*e_X) + cross( (betadot*e_pin + psidot*e_X), r_CD );
RHS = -vD*e_Z;

% Assemble and solve the velocity equation
velocityEquation = LHS - RHS;
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), velocityEquation(3), ...
    vD, betadot, psidot );

vD_computed = double( velocitySolutionStruct.vD );
betadot_computed = double( velocitySolutionStruct.betadot );
psidot_computed = double( velocitySolutionStruct.psidot );

% Print results
fprintf( 'v_D      = %0.4f m/s\n', vD_computed );
fprintf( 'thetadot = %0.4f rad/s\n', betadot_computed );
fprintf( 'psidot   = %0.4f rad/s\n\n', psidot_computed );
    