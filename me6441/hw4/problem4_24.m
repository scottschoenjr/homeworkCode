% Problem 4-24 in _Engineering Dynamics_ (Ginsberg, 2008)

clear all
close all
clc

% Velocity

% Declare symbolic variables
syms vA betadot psidot

% Define the vector connecting the ball and socket to the collar
gamma  = 75.*pi./180; % Angle that A's guide makes with y-axis
yA = 3.*(1 - cos(gamma) );
zA = 3.*sin(gamma);
rA = [ 0; yA; zA ]; % Collar A [m]
rB = [ 2; 0; 0 ];   % Collar B [m]
r_AB = rA - rB;     % Difference

% A bit of geometry to find the current angle 
Gamma = sqrt( 3.^(2) + 2.^(2) );
L = norm(r_AB);
psi = acos( Gamma./L );

% Define the unit vectors
e_X = [1; 0; 0];
e_Y = [0; 1; 0];
e_Z = [0; 0; 1];

% Define the unit vector for the pin
e_pin = cross( e_X, r_AB );
e_pin = e_pin./norm(e_pin);

% Define the velocity along the x axis of collar B
vB = 30*e_X;

% Set up velocity equation to solve
LHS = vB + cross( (betadot*e_pin + psidot*e_X), r_AB );
RHS = -vA*cos(gamma)*e_Y + vA*sin(gamma)*e_Z;
velocityEquation = LHS + RHS;

% Run the solver
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), velocityEquation(3), ...
    vA, betadot, psidot );

vA_computed = double( velocitySolutionStruct.vA );
betadot_computed = double( velocitySolutionStruct.betadot );
psidot_computed = double( velocitySolutionStruct.psidot );

% Print results
fprintf( 'v_A     = %0.4f m/s\n', vA_computed );
fprintf( 'detadot = %0.4f rad/s\n', betadot_computed );
fprintf( 'psidot  = %0.4f rad/s\n\n', psidot_computed );

% Accleration
syms aA betaddot psiddot

% Known accleration of B
aB = -500*e_X;

% Define rotational velocity we just solved for
omegaAB = betadot_computed*e_pin + psidot_computed*e_X;

% Define angular acceleration
alpha = ...
    betaddot*e_pin + betadot_computed*cross( omegaAB, e_pin ) + psiddot*e_X;

% Assemble the acceleration equation
LHS = aB + cross( alpha, r_AB ) + cross( omegaAB, cross( omegaAB, r_AB ) );
RHS = aA*cos(gamma)*e_Y - aA*sin(gamma)*e_Z;
accelerationEquation = LHS - RHS;

accelerationSolutionStruct = solve( ...
    accelerationEquation(1), accelerationEquation(2), accelerationEquation(3), ...
    aA, betaddot, psiddot ... % Variables to solve for
    );

aA = double( accelerationSolutionStruct.aA );
betaddot = double( accelerationSolutionStruct.betaddot );
psiddot = double( accelerationSolutionStruct.psiddot );

% Print results
fprintf( 'a_A      = %0.4f m/s^2\n', aA );
fprintf( 'betaddot = %0.4f rad/s^2\n', betaddot );
fprintf( 'psiddot  = %0.4f rad/s^2\n', psiddot );


    