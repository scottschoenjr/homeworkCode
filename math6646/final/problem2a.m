% Solve for coefficients

clear all
close all
clc

% Set up matrix
A = [ ...
       1, 1,  0,    0; ...
       1, 0, -1,   -1; ...
    1./2, 0, -1,    0; ...
    1./6, 0, -1./2, 0 ...
    ];
b = [-1; -2; -2; -4./3];
A\b