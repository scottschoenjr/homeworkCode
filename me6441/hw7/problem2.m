% Problem 2 of HW 7

clear all;
close all;
clc; 

% Define plate parameters
m = 0.2; % [kg]
l_x = 90E-3; % [m]
l_y = 60E-3; % [m]
L = 100E-3; % [m]

% Define angle
theta = atan2( l_y./2, l_x./2 );

% Define unit vectors
e_x1 = [1; 0; 0];
e_y1 = [0; 1; 0];
e_z1 = [0; 0; 1];

e_x = cos(theta).*e_x1 - sin(theta).*e_y1;
e_y = sin(theta).*e_x1 + cos(theta).*e_y1;
e_z = e_z1;

% Define angular velocities
Omega = 100; % [rad/s]
omega = Omega.*e_x1;
alpha = 0.*omega;

% Define applied torque
Gamma0 = 3E3; % [Nm]

% Define omega matrix
omega_x = dot( omega, e_x );
omega_y = dot( omega, e_y );
omega_z = dot( omega, e_z );
omegaVec = [omega_x; omega_y; omega_z];
omegaMat = [ ...
           0, -omega_z,  omega_y; ...
     omega_z,        0, -omega_x; ...
    -omega_y,  omega_x,        0 ...
    ];


% Define inertia matrix parameters
I_xx = m.*l_x.^(2)./12;
I_yy = m.*l_y.^(2)./12;
I_zz = m.*(l_x.^(2) + l_y.^(2))./12;
I = [ I_xx,    0,   0; ...
         0, I_yy,   0; ...
         0,    0, I_zz ...
     ];
 
% Define moment vector
M = I*alpha + omegaMat*I*omegaVec;

% Define total torque in z
A_z = m*Omega.^(2)*cos(theta).*sin(theta).*(l_y.^(2) - l_x.^(2))./(24.*L)

% Define angular acceleration components
alpha_x = 24.*(A_z.*sin(theta) - Gamma0.*cos(theta))./(m.*l_x.^(2))
alpha_y = 24.*(A_z.*cos(theta) + Gamma0.*sin(theta))./(m.*l_y.^(2)) 
alpha = alpha_x.*cos(theta) + alpha_y.*sin(theta)


