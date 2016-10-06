% Matlab starter code
% Dr. Ferri

L = 1; %m
g = 10; %m/s^2

% create a 2x1 "anonymous function" to describe the 
% right-hand-side of the state equations, f(t,x)

fnc = @(t,x)  [x(2); -(g/L)*sin(x(1))];

t0 = 0; tf = 20; % sec
x0 = [1; -0.5]; % initial conditions of 1 rad and -0.5 rad/s

% set "options" to prescribe tolerances for accuracy.
% Reduce reltol, abstol, and/or maxstep  for more accurate results
options = odeset('reltol',1e-7,'abstol',1e-7,'maxstep',(tf-t0)/100);

[T,X]= ode45(fnc,[t0 tf],x0,options);

% T is a px1 column vector of times (ranging from t0 to tf)
% X is a px2 matrix; the first column gives x1, the second x2, etc

theta = X(:,1);  % extract first column
thetadot = X(:,2);  % extract 2nd column


figure(1)
plot(T,theta)
xlabel('time (s)')
ylabel('angle (rad)')

figure(2)
plot(T,thetadot)
xlabel('time (s)')
ylabel('angular velocity (rad/s)')