% Problem 4-14 in _Engineering Dynamics_ (Ginsberg, 2008)

clear all
close all
clc

% Define symbolic values to solve for
syms vAB omegaBC

% Compute values of rAB and theta for the known beta (side lengths are
% normalized to L)
L = 1;
armLength = 0.8; % Relative to L
beta = 20.*(pi./180);
rAB = L.^(2) + armLength.^(2) - 2.*L.*armLength.*cos( beta ); % LOC
theta = atan( 0.8*sin(beta)./(1 - 0.8*cos(beta)) ); % LOS

% We'll normalize all rotational velocities to omegaA 
omegaA = 1;

% Define the inertial unit vectors 
e_X = [1; 0; 0];
e_Y = [0; 1; 0];
e_Z = [0; 0; 1];

% Define the radial unit nomal vector
e_r = cos(theta)*e_X + sin(theta)*e_Y;

% Define BC quantities
rBC = -armLength.*cos(beta)*e_X + armLength.*sin(beta)*e_Y;
betadot = -omegaBC*e_Z;

% Define the velocity of the point in terms of the disc
LHS = vAB*e_r ...
    + omegaA*rAB*cos(theta)*e_Y - omegaA*rAB*sin(theta)*e_X;
RHS = cross( betadot, rBC );

% Solve the equation when it's 0 for omegaBC and vAB
velocityEquation = LHS - RHS;
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), ...
    vAB, omegaBC );

vAB_computed = double( velocitySolutionStruct.vAB );
omegaBC_computed = double( velocitySolutionStruct.omegaBC );

% Print results
fprintf( 'v_AB     = %0.4f m/s\n', vAB_computed );
fprintf( 'thetadot = %0.4f *omega_A rad/s\n', omegaBC_computed );

% Acceleration

syms a_B alpha_BC

% Define acceleration terms
vBrel = vAB_computed*e_r;
aBrel = -rAB*omegaA.^(2)*e_r;

% Define left and right-hand sides
omegaBC = omegaBC_computed*e_Z;
omegaAB = omegaA*e_Z;
rABvec = rAB*e_r;
LHS = cross( alpha_BC*e_Z, rBC ) + cross( omegaBC, cross( omegaBC, rBC ) );
RHS = vAB_computed*e_r + cross( omegaAB, cross( omegaAB, rABvec) );


% Solve
accelerationEquation = RHS - LHS;
accelerationSolutionStruct = solve( ...
    accelerationEquation(1), accelerationEquation(2), ...
    a_B, alpha_BC );

aB_computed = double( accelerationSolutionStruct.a_B );
alphaBC_computed = double( accelerationSolutionStruct.alpha_BC );

% Print results
fprintf( 'a_B     = %0.4f m/s^2\n', aB_computed );
fprintf( 'alphaBC = %0.4f *omega^2 rad/s^2\n', alphaBC_computed );





    