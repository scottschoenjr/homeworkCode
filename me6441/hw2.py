# Code for Homework 2, ME-6441

import numpy as np

# Useful constants
deg2rad = np.pi/180;
rad2deg = 1/deg2rad;

# Define rotation matrices
def RotX(theta):
    R = np.matrix([[1,              0,            0 ], 
                   [0,  np.cos(theta), np.sin(theta)], 
                   [0, -np.sin(theta), np.cos(theta)]])
    return R
#
def RotY(theta):
    R = np.matrix([[np.cos(theta), 0, -np.sin(theta)], 
                   [0,             1,            0  ], 
                   [np.sin(theta), 0,  np.cos(theta)]]) 
    return R
#
def RotZ(theta):
    R = np.matrix([[ np.cos(theta), np.sin(theta), 0],  
                   [-np.sin(theta), np.cos(theta), 0],
                   [             0,             0, 1]]) 
    return R
#

# Problem 1 -----------------------------------------------------------------
r = np.matrix([[30E-3], [80E-3], [50E-3]])
Ry = RotY(65*deg2rad);
Rx = RotX(-135*deg2rad);
R = np.dot(Rx, Ry); # Total rotation matrix

# Part a
r1 = np.dot(R, r);
print("Problem 1(a): r = " )
print( r1 )
print("")

# Part b
r2 = np.dot(R.T, r);
print("Problem 1(b): r = " )
print(r2)
print("")

# Problem 2 -------------------------------------------------------------------
d = 120E-3;

L0 = 250E-3;
theta0 = 0;
beta0 = np.pi/2;
gamma0 = 0;

L1 = 500E-3;
theta1 = 2*np.pi/3;
beta1 = 2*np.pi/3;
gamma1 = -np.pi/2;

# First position
mat1 = RotZ(theta0);
mat2 = RotY(beta0);
mat3 = RotX(gamma0);
Rtot = np.dot( mat3, np.dot(mat2, mat1) );
r = np.matrix([[L0],[0],[d]])
x0 = np.dot(Rtot, r)

# Second position
mat1 = RotZ(theta1);
mat2 = RotY(beta1);
mat3 = RotX(gamma1);
Rtot = np.dot( mat3, np.dot(mat2, mat1) );
r = np.matrix([[L1],[0],[d]])
x1 = np.dot(Rtot, r)

# Print results
print("Problem 2: Delta r = " )
print( x1 - x0 )
print("" )

# Problem 4 -------------------------------------------------------------------

# Angles and distances of interest
theta = 120*deg2rad;
phi = -70*deg2rad;
AB = 400E-3;
BC = 400E-3;
CD = 200E-3;

# Define rotation matrices
R0i = RotX(theta);
R1i = np.dot( RotZ(50*deg2rad), R0i );
R2i = np.dot( RotX(phi), R1i );

# De
v0 = np.matrix([[AB], [0], [ 0]]);
v1 = np.matrix([[BC], [0], [ 0]]);
v2 = np.matrix([[ 0], [0], [CD]]);

R = np.dot( R0i.T, v0 ) + np.dot( R1i.T, v1 ) + np.dot( R2i.T, v2 );

# Print results
print("Problem 4: R = " )
print( R )
print("" )

