% Problem 4-10 in "Engineering Dynamics" (Ginsberg, 2008)

% Define symbolic values to solve for
% syms vB betadot

% Compute values of L and D (lengths are normalized to R)
R = 1;
D = 1.75.*R;
L = 1.25.*R;

% Compute theta at this instant
theta = asin( R./L );

% We'll normalize all rotational velocities to omega (=thetadot) 
omega = 1;

% Define the inertial unit vectors 
e_X = [1; 0];
e_Y = [0; 1];

% Define the velocity of the collar in terms of arm AB
LHS = omega*(D - L*cos(theta))*e_Y - omega*L*sin(theta)*e_X;

RHS = - omegaBC*L*(1 - armLength*cos(beta))*e_Y ...
    + omegaBC*armLength*L*sin(beta)*e_X;

% Solve the equation when it's 0 for omegaBC and vAB
velocityEquation = LHS - RHS;
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), ...
    vAB, omegaBC );

vAB_computed = double( velocitySolutionStruct.vAB );
omegaBC_computed = double( velocitySolutionStruct.omegaBC );

% Print results
fprintf( 'v_AB      = %0.4f m/s\n', vAB_computed );
fprintf( 'thetadot = %0.4f rad/s\n', omegaBC_computed );
    