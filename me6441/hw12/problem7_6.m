% Ex 7.6

clear all
close all
clc

% -------------------------
% Generalized coordinates
% q1  - psi
% dq1 - psiDot
% q2  - theta
% dq2 - thetaDot
% -------------------------

% Define symbolic values
syms q1 dq1 q2 dq2 kT I1 L m Gamma P

% Define unit vectors
eX = [ 1; 0; 0 ];
eY = [ 0; 1; 0 ];
eZ = [ 0; 0; 1 ];

% Kinetic Energy -----

% Of gimbal
HG = I1.*dq1.*eX;
omega = dq1.*eX;
T1 = (1./2).*omega.*HG;

% Of arm
IXX = 2.*m.*L.^(2).*(sin(q2)).^(2);
IZZ = 2.*m.*L.^(2);
HG = IXX.*dq1.*eX + IZZ.*dq2.*eZ;
omega = dq1.*eX + dq2.*eZ;
T2 = (1./2).*omega.*HG;

% Total KE
T = T1 + T2;

% Potential Energy -----

% Spring
V = (1./2).*kT.*q2.^(2);

% Get equation of motion with fulldiff
eomPsi = fulldiff( diff(T, dq1), {q1, q2}) ...
    - diff(T, q1) + diff(V, q1)
eomTheta = fulldiff( diff(T, dq2), {q1, q2}) ...
    - diff(T, q2) + diff(V, q2)
