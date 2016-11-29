% Ex 7.30
% Mass is given as a mass-per-unit-length, sigma.
%    mass of shorter bar is sL = sigma*L;
%    mass of longer bar is m2 = sigma*(3*L/2) = 3*sL/2

syms q1 dq1 k L m g P Q1 sL;

% Kinetic Energy
IA = (sL*L^2)/3;  % pure rotation about A
L2 = 3*L/2; m2 = sL*3/2;
IG = (m2*L2^2)/12;
rGA = [5*L*cos(q1)/4; 3*L*sin(q1)/4; 0];
vG = diff(rGA,q1)*dq1;
T1 = (IA*dq1^2)/2;
T2 = (m2*vG.'*vG)/2 + (IG*dq1^2)/2;  % Note: a.'*b = transpose(a)*b
T = T1+T2;

% Potential Energy
Vg = sL*g*L*sin(q1)/2 + m2*3*L*sin(q1)/4;
DEL = 2*L*cos(q1) - 2*L*cos(pi/4);
Vs = (k*DEL^2)/2;
V = Vg + Vs;

% Virtual work
rDA = [L*cos(q1)/2; 3*L*sin(q1)/2; 0];
del_rD = diff(rDA,q1);
Pvec = P*[0; -1; 0];  % vector expression for P in downward (-j) direction
Q1 = Pvec.'*del_rD;  % dot product of P with del_rD

% All time varying variables go where {q1} is.  Note that if 'q1' is
% specified as an argument for fulldiff, then dq1 is assumed to be time
% dependent as well.
eq_q1 = fulldiff(diff(T,dq1),{q1}) - diff(T,q1) + diff(V,q1) - Q1;

% EOM can be written:
%
%  d2q1*(sL*L^2)*(35/24  + (3/2)*sin(q1)^2) 
%      + (3/2*(sL*L^2)*sin(q1)*cos(q1)*dq1^2
%      + 2*k*(L^2)*(sqrt(2) - 2*cos(q1))*sin(q1)
%      + (13/8)*sL*L*g*cos(q1) = -(3/2)*P*L*cos(q1)