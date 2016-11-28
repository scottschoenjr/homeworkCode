% Ex 7.30

clear all
close all
clc

% -------------------------
% Generalized coordinates
% q1  - theta
% dq1 - thetaDot
% -------------------------

% Define symbolic values
syms q1 dq1 k L m g P Q1 sL

% Define unit vectors
eX = [ 1; 0; 0 ];
eY = [ 0; 1; 0 ];
eZ = [ 0; 0; 1 ];

% Kinetic Energy -----

% Arm 1
m1 = m;
IA = (m1.*L^(2))/3;
T1 = ( IA.*dq1.^(2) )./2;

% Arm 2
L2 = 3.*L/2; 
m2 = 3/2.*m1;
IG = (m2.*L2^2)/12; % About center of mass
rGA = (5./4)*L*cos(q1).*eX +  (3/4).*L.*sin(q1).*eY;
vG = diff(rGA, q1).*dq1; % Velocity of center of mass
T2 = (1./2).*m2*vG.'*vG + (1./2).*IG.*dq1.^(2);

% Total KE
T = T1 + T2;

% Potential Energy -----

% Gravitational
Vg = m1.*g.*(L./2).*sin(q1) + m2.*g.*(3.*L/4).*sin(q1);

% Spring
deltaX = 2.*L.*cos(q1) - 2*L*cos(pi./4);
Vs = (1./2).*k.*deltaX.^(2);

% Total PE
V = Vg + Vs;

% Virtual work
rDA = (L/2).*cos(q1).*eX + (3.*L./2)*sin(q1).*eY;
delta_rD = diff(rDA, q1);
Pvec = -P.*eY; % Applied force
Q1 = Pvec.'*delta_rD; % Generalized force

% Get equation of motion with fulldiff
eom = fulldiff( diff(T, dq1), {q1}) ...
    - diff(T,q1) + diff(V, q1) - Q1

% EOM can be written:
%
%  d2q1*(sL*L^2)*(35/24  + (3/2)*sin(q1)^2) 
%      + (3/2*(sL*L^2)*sin(q1)*cos(q1)*dq1^2
%      + 2*k*(L^2)*(sqrt(2) - 2*cos(q1))*sin(q1)
%      + (13/8)*sL*L*g*cos(q1) = -(3/2)*P*L*cos(q1)
