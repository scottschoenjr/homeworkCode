% Problem 4-23 in "Engineering Dynamics" (Ginsberg, 2008)

clear all
close all
clc

% Velocity

% Declare symbolic variables
syms vD omegaCD

% Define the vector connecting the ball and socket to the collar
rC = [ 0; 3; 0 ];         % Cap C [m]
rD = [ -0.5; 0; 2.5981 ]; % Collar D [m]
r_CD = rC - rD;

% Define the unit vectors
e_x = [1; 0; 0];
e_y = [0; 1; 0];
e_z = [0; 0; 1];

% Define the unit vector for the pin
e_pin = cross( e_Y, r_CD );
e_pin = e_pin./norm(e_pin);
    