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
    return R;
#
def RotY(theta):
    R = np.matrix([[np.cos(theta), 0, -np.sin(theta)], 
                   [0,             1,            0  ], 
                   [np.sin(theta), 0,  np.cos(theta)]]) 
    return R;
#
def RotZ(theta):
    R = np.matrix([[ np.cos(theta), np.sin(theta), 0],  
                   [-np.sin(theta), np.cos(theta), 0],
                   [             0,             0, 1]]) 
    return R;
#

# Problem 1 -----------------------------------------------------------------
r = np.matrix([[30], [80], [50]])
Ry = RotY(65*deg2rad);
Rz = RotZ(-145*deg2rad);
R = np.dot(Rz, Ry); # Total rotation matrix

# Part a
r1 = np.dot(R, r);
print("Problem 1(a) r =  ");
print( r1 );
print( "" );

# Part b
r2 = np.dot(R.T, r);
print("Problem 1(b) r = " );
print(r2);
print("");

# Problem 2 -------------------------------------------------------------------
d = 120;

L0 = 250;
theta0 = 0;
beta0 = np.pi/2;
gamma0 = 0;

L1 = 500;
theta1 = 2*np.pi/3;
beta1 = 2*np.pi/3;
gamma1 = -np.pi/2;

# First position
mat1 = RotZ(theta0);
mat2 = RotX(np.pi/2 - beta0);
mat3 = RotY(gamma0);
Rtot = np.dot( mat3, np.dot(mat2, mat1) );
r = np.matrix([[0],[L0],[d]])
r0 = np.dot(Rtot.T, r)

# Second position
mat1 = RotZ(theta1);
mat2 = RotX(np.pi/2 - beta1);
mat3 = RotY(gamma1);
Rtot = np.dot( mat3, np.dot(mat2, mat1) );
r = np.matrix([[0],[L1],[d]])
r1 = np.dot(Rtot.T, r)

# Print results
print("Problem 2 Delta r = " )
print( r1 - r0 )
print( "" )

# Problem 4 -------------------------------------------------------------------

# Angles and distances of interest
theta = 120*deg2rad;
phi = -70*deg2rad;
AB = 400;
BC = 400;
CD = 200;

# Define rotation matrices
R1i = RotX(theta);
R2i = np.dot( RotZ(np.pi - 50*deg2rad), R1i );
R3i = np.dot( RotX(-phi), R2i );

# Define the position vectors in each frame
v1 = np.matrix([[AB], [0],   [0]]);
v2 = np.matrix([[BC], [0],   [0]]);
v3 = np.matrix([[ 0], [-CD], [0]]);

# Compute the vector in the inertial frame
R = np.dot( R1i.T, v1 ) + np.dot( R2i.T, v2 ) + np.dot( R3i.T, v3 );

# Print results
print("Problem 4 R =  ")
print( R )
print("" )