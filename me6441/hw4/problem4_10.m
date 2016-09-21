% Problem 4-10 in _Engineering Dynamics_ (Ginsberg, 2008)

clear all
close all
clc

% Define symbolic values to solve for
syms u betadot

% Compute values of L and D (lengths are normalized to R)
R = 1;
D = 1.75.*R;
L = 1.25.*R;

% Compute theta at this instant
theta = asin( R./L );

% We'll normalize all rotational velocities to omegaAB (= thetadot) 
omegaAB = 1;

% Define the inertial unit vectors 
e_X = [1; 0];
e_Y = [0; 1];

% Define the velocity of the collar in terms of arm AB
LHS = 0.75*omegaAB*R*e_Y + omegaAB*R*e_X;
RHS = (u - R*betadot)*e_X + R*betadot*e_Y;

% Solve the equation when it's 0 for omegaBC and vAB
velocityEquation = LHS - RHS;
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), ...
    u, betadot );

u_computed = double( velocitySolutionStruct.u );
betadot_computed = double( velocitySolutionStruct.betadot );

% Print results
fprintf( 'u       = %0.4f m/s\n', u_computed );
fprintf( 'betadot = %0.4f *omega rad/s\n\n', betadot_computed );

% Acceleration
syms udot betaddot

% Define LHS and RHS of accleration expression
LHS = omegaAB.^(2)*L*( cos(theta)*e_X - sin(theta)*e_Y );
RHS = udot*e_X - (u_computed.^(2)/R)*e_Y + R*betaddot*( e_X - e_Y ) - ...
    R*betadot_computed.^(2)*(e_X + e_Y) + 2*betadot_computed*u_computed*e_Y;
accelerationEquation = RHS - LHS;

% Solve
accelerationSolutionStruct = solve( ...
    accelerationEquation(1), accelerationEquation(2), ...
    udot, betaddot );

udot_computed = double( accelerationSolutionStruct.udot );
betaddot_computed = double( accelerationSolutionStruct.betaddot );

% Print results
fprintf( 'udot     = %0.4f m/s^2\n', udot_computed );
fprintf( 'betaddot = %0.4f *omega rad/s^2\n', betaddot_computed );

    