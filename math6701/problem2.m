% Problem 2

clear all
close all
clc

A = [ ...
    2, 2, 0; ...
    4, 0, 0; ...
    1, 2, 1 ...
    ];

% Find eigenvalues and vectors
[V, D] = eig(A);

lambda1 = D(1, 1);
lambda2 = D(2, 2);
lambda3 = D(3, 3);  

v1 = V(:, 1);
v2 = V(:, 2);
v3 = V(:, 3);

% Coupute the normal way
A10_1 = V*(D^10)*inv(V);

% Compute coefficients
m = [ ...
    1, lambda1, lambda1.^(2); ...
    1, lambda2, lambda2.^(2); ...
    1, lambda3, lambda3.^(2)  ...
    ];

n =10;
L = [ ...
    lambda1.^(n); lambda2.^(n); lambda3.^(n) ...
    ];

% Solve for c0
c = m\L;

% Compute A10 the second way
A10_2 = c(1).*eye(3) + c(2).*A + c(3).*A^(2);
    

