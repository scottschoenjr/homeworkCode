% ME 6441 HW 12 Problem 4 

clear all
close all
clc

% -------------------------
% Generalized coordinates
% q1  - psi
% dq1 - psiDot
% q2  - theta
% dq2 - thetaDot
% q3  - q
% dq3 - qDot
% -------------------------

% Define symbolic values
syms q1 dq1 q2 dq2 q3 dq3 L m1 m2 F M g k

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

% Rotational velocity of rod
Izz = (1./12).*m1.*L.^(2);
Iyy = (1./12).*m1.*L.^(2);
HG = -Iyy.*dq2.*ey + Izz.*dq1.*ez;
omega = -dq2.*ey + dq1.*eZ;
T2 = (1./2).*omega.'*HG;

% Translational velocity of collar
rC = q3.*ex;
vC = dq3.'*ex + cross( omega, rC );
T1 = (1./2).*m2.*vC.'*vC;

% Total KE
T = T1 + T2;
simplify(T)

% Potential Energy -----
V_AB = -m1*g*( (L./2).*sin(q2) ); % Rod (just gravitational)
V_C = -m2*g*( q3.*sin(q2) ) + (1./2).*k.*(q3 - (L./2)); % Collar
V = V_AB + V_C;

% Get equations of motion with fulldiff
eomTheta = fulldiff( diff(T, dq2), {q2, q3}) ...
    - diff(T, q2) + diff(V, q2);
pretty( simplify( eomTheta ) )
eomQ = fulldiff( diff(T, dq3), {q2, q3}) ...
    - diff(T, q3) + diff(V, q3);
pretty( simplify( eomQ ) )