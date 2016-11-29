% ME 6441 HW 12 Problem 3 

clear all
close all
clc

% -------------------------
% Generalized coordinates
% q1  - psi
% dq1 - psiDot
% q2  - theta
% dq2 - thetaDot
% q3  - xi
% dq3 - xiDot
% -------------------------

% Define symbolic values
syms q1 dq1 q2 dq2 q3 dq3 L m F g

% Define inertial unit vectors
eX = [ 1; 0; 0 ];
eY = [ 0; 1; 0 ];
eZ = [ 0; 0; 1 ];

% Define unit vectors attached to rod
gamma = pi/2 - q2;
Ry = [ cos(gamma), 0, sin(gamma); 
                0, 1,          0;
      -sin(gamma), 0, cos(gamma) ];
ex = Ry*eX;
ey = Ry*eY;
ez = Ry*eZ;

% Kinetic Energy -----

% Translational velocity
rG = (q3 - L/2).*ex;
vG = diff( rG, q1 )*dq1 + diff( rG, q2 )*dq2 + diff( rG, q3 )*dq3;
T1 = (1./2).*m.*vG.'*vG;

% Rotational velocity
Izz = (1./3).*m.*L.^(2);
Iyy = (1./3).*m.*L.^(2);
omega = -dq2.*ey + dq1.*eZ;
omega_z = omega.'*ez;
omega_y = omega.'*ey;
HG = -Iyy.*omega_y.*ey + Izz.*omega_z.*ez;
T2 = (1./2).*omega.'*HG;

% Total KE
T = T1 + T2;
simplify(T)

% Potential Energy -----
V = m*g*( rG.'*eZ );

% Get equations of motion with fulldiff
eomPsi = fulldiff( diff(T, dq1), {q1, q2, q3}) ...
    - diff(T, q1) + diff(V, q1)
eomTheta = fulldiff( diff(T, dq2), {q1, q2, q3}) ...
    - diff(T, q2) + diff(V, q2)
eomXi = fulldiff( diff(T, dq1), {q1, q2, q3}) ...
    - diff(T, q2) + diff(V, q2)
