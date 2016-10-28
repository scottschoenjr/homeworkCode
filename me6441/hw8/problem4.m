% ME 6441 Homework 8 Problem 4

clear all
close all
clc

P = 16; % [N]
kc = 2; % [m]
m = 2; % [kg]
Ro = 3; % [m]
Ri = 2; % [m]

% Part a
alpha = ( P.*(Ri - Ro) )./( m.*(Ri.^(2) + kc.^(2)) );
a = Ri*alpha

% Part b
aCx = P/m
alpha = -P*Ro/(m*kc.^(2))