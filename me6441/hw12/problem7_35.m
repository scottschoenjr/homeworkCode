% Exercise 7.35, 3rd edition

%-------------------------
% (Mandatory) Notation
%-------------------------
% q1  = x
% dq1
% q2  = theta
% dq2
%-------------------------

syms q1 dq1 q2 dq2 m1 m2 g L k ;

I = [1; 0; 0]; J = [0; 1; 0]; K = [0; 0; 1];

% Kinetic energy 
IG = (m2*L^2)/12;   
rGO = (L/2)*[sin(q2); -cos(q2); 0];
vG = dq1*I + cross(dq2*K,rGO) ;
T = (m1*dq1^2)/2 + (m2*vG.'*vG)/2 + (IG*dq2^2)/2;
T = simplify(T)

% T =
%  
% (dq1^2*m1)/2 + (dq1^2*m2)/2 + (L^2*dq2^2*m2)/6 + (L*dq1*dq2*m2*cos(q2))/2

% Potential Energy
V = (2*k*q1^2)/2 + -m2*g*(L/2)*cos(q2)

% Virtual work of vertical (upward) force applied to end B

Q1 = 0;
Q2 = 0;

% All time varying variables go where {q1} is.  Note that if 'q1' is
% specified as an argument for fulldiff, then dq1 is assumed to be time
% dependent as well.

disp('eq1')
eq_q1 = fulldiff(diff(T,dq1),{q1,q2}) - diff(T,q1) + diff(V,q1) - Q1;
eq_q1 = simplify(eq_q1)

disp(' ')
disp('eq2')
eq_q2 = fulldiff(diff(T,dq2),{q1,q2}) - diff(T,q2) + diff(V,q2) - Q2;
eq_q2 = simplify(eq_q2)

% eq_q1 =
%  
% 2*k*q1 + d2q1*(m1 + m2) + (L*d2q2*m2*cos(q2))/2 
%           - (L*dq2^2*m2*sin(q2))/2 = 0
%  
% eq_q2 =
%  
% (L*m2*(2*L*d2q2 + 3*d2q1*cos(q2) + 3*g*sin(q2)))/6 = 0