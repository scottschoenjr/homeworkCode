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
e_X = [1; 0];
e_Y = [0; 1];

% Define the velocity of the point in terms of the disc
LHS = vAB*( cos(theta)*e_X + sin(theta)*e_Y ) ...
    + omegaA*rAB*cos(theta)*e_Y - omegaA*rAB*sin(theta)*e_X;
RHS = omegaBC*armLength*L*cos(beta)*e_Y ...
    + omegaBC*armLength*L*sin(beta)*e_X;

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
    