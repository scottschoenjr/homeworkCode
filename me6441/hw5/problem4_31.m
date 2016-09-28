% Problem 4-31 in _Engineering Dynamics_ (Ginsberg, 2008)

clear all
close all
clc

% Define symbolic values to solve for
syms omegaC omegaAB

% Set values of collar speed and disc radius (results will be normalized to
% these)
R = 1;
v = 1;

% Define the inertial unit vectors 
e_X = [1; 0];
e_Y = [0; 1];

% Define the velocity of the collar in terms of arm AB
LHS = (v - 1.5*R*omegaAB)*e_X - 2*R*omegaAB*e_Y;
RHS = -R*omegaC*e_X + 0.5*R*omegaC*e_Y;

% Solve the equation when it's 0 for omegaBC and vAB
velocityEquation = LHS - RHS;
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), ...
    omegaAB, omegaC );

omegaAB_computed = double( velocitySolutionStruct.omegaAB );
omegaC_computed = double( velocitySolutionStruct.omegaC );

% Print results
fprintf( 'omega_AB = %0.4f*(v/R)\n', omegaAB_computed );
fprintf( 'omega_C  = %0.4f*(v/R)\n\n', omegaC_computed );

% Acceleration
syms omegadotAB omegadotC

% Define left- and right-hand sides
LHS = ( 1.5*R*omegadotAB + 2*omegaAB_computed^(2) )*e_X + ...
    ( 2*R*omegadotAB - 1.5*R*omegaAB_computed^(2) )*e_Y;
RHS = ( -R*omegadotC - 0.5*omegaC_computed^(2) )*e_X + ...
    ( 0.5*R*omegadotC )*e_Y;

% Solve acceleration equation
accelerationEquation = LHS + RHS;
accelerationSolutionStruct = solve( ...
    accelerationEquation(1), accelerationEquation(2), ...
    omegadotAB, omegadotC );

omegadotAB_computed = double( accelerationSolutionStruct.omegadotAB );
omegadotC_computed = double( accelerationSolutionStruct.omegadotC );

% Print results
fprintf( 'omegadot_AB = %0.4f*(v/R)^2\n', omegadotAB_computed );
fprintf( 'omegadot_C  = %0.4f*(v/R)^2\n\n', omegadotC_computed );

