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
omegaA = 2.*pi.*(900./60); % [rad/s]
rBC = 1; % [m]

% Define the unit vector for the pin
e_pin = cross( e_Y, r_CD );
e_pin = e_pin./norm(e_pin);

% Define the sides of the velocity equation
omegaCD = (betadot*e_pin + psidot*e_Y + omegaA*e_Z);
LHS = -(rBC*omegaA*e_X) + cross( omegaCD, r_CD );
RHS = vD*e_Z;

% Assemble and solve the velocity equation
velocityEquation = LHS - RHS;
velocitySolutionStruct = solve( ...
    velocityEquation(1), velocityEquation(2), velocityEquation(3), ...
    vD, betadot, psidot );

vD_computed = double( velocitySolutionStruct.vD );
betadot_computed = double( velocitySolutionStruct.betadot );
psidot_computed = double( velocitySolutionStruct.psidot );
omegaCD_computed = betadot_computed*e_pin + psidot_computed*e_Y + omegaA*e_Z;

% Print results
fprintf( 'v_D     = %0.4f m/s\n', vD_computed );
fprintf( 'betadot = %0.4f rad/s\n', betadot_computed );
fprintf( 'psidot  = %0.4f rad/s\n', psidot_computed );
fprintf( 'omegaCD = [%0.4f, %0.4f, %0.4f] \n\n', omegaCD_computed );

% Acceleration

% Declare variables
syms aD betaddot psiddot

% Define the sides of the equation
omegaCD = betadot_computed*e_pin + psidot_computed*e_Y + omegaA*e_Z;
alpha = betaddot*e_pin + betadot_computed*cross( omegaCD, e_pin ) + ...
    psiddot*e_Y - psidot_computed*cross( omegaA*e_Z, e_Y );
RHS = rBC.*omegaA.^(2)*e_Y + ...
    cross( alpha, r_CD ) + ...
    cross( omegaCD, cross( omegaCD, r_CD ) );
LHS = aD*e_Z;

% Define and solve accleration equation
accelerationEquation = RHS - LHS;
accelerationSolutionStruct = solve( ...
    accelerationEquation(1), accelerationEquation(2), accelerationEquation(3), ...
    aD, betaddot, psiddot ...
    );

aD = double( accelerationSolutionStruct.aD );
betaddot = double( accelerationSolutionStruct.betaddot );
psiddot = double( accelerationSolutionStruct.psiddot );
alphaC_computed = betaddot*e_pin + betadot_computed*cross( omegaCD, e_pin ) + psiddot*e_Y;

% Print results
fprintf( 'a_D       = %0.4f m/s^2\n', aD );
fprintf( 'betaddot  = %0.4f rad/s^2\n', betaddot );
fprintf( 'psiddot   = %0.4f rad/s^2\n', psiddot );
fprintf( 'alphaC = [%0.4f, %0.4f, %0.4f] \n\n', alphaC_computed );
    