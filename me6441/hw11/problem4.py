# -*- coding: utf-8 -*-
"""
Problem 4 Homework 11 ME6441
Exercize 7.36 in Ginsberg's Advanced Dynamics
"""

import sympy
import sympy.physics.mechanics as mech

# Set dynamics variables
phi = mech.dynamicsymbols('theta'); # Blade angle
phiDot = mech.dynamicsymbols('thetaDot'); # Blade angluar velocity

# Create symboloic variables we'll need
m, L, kT, Omega, epsilon, Izz = sympy.symbols('m L kT Omega epsilon Izz');

# Create inertial frame
inertialFrame = mech.ReferenceFrame('inertialFrame');

# Create reference frame at A that rotates with Omega
N = mech.ReferenceFrame( 'N' );
N.set_ang_vel(inertialFrame, -Omega*inertialFrame.z);

# Set point A
A = mech.Point('A');
A.set_vel(N, 0);

# Set point B
B = mech.Point( 'B' );
B.set_vel( N, 0 ); # B has no velocity in the rotating frame N
B.set_pos( A, epsilon*N.x ); # phi measured from x-axis along AB

# Create frame for rod BC 
rodBCFrame = mech.ReferenceFrame('rodBCFrame');
rodBCFrame.set_ang_vel( N, -phiDot*N.z );

# Create point for the center of mass of BC
comBC_x = (L/2)*sympy.cos(phi);
comBC_y = (L/2)*sympy.sin(phi);
comBC_vx = phiDot*comBC_x;
comBC_vy = phiDot*comBC_y;

comBC = mech.Point('comBC');
comBC.set_vel( N, comBC_vx*N.x + comBC_vy*N.y ); 
comBC.set_pos( B, comBC_x*N.x + comBC_y*N.y );

# Create rigid body for blade BC
I_BC = mech.inertia(N, 0, 0, Izz);
rodBC = mech.RigidBody( 'rodBC', comBC, rodBCFrame, m, (I_BC, comBC) );

# Get the kinetic energy of each component to check results
mech.kinetic_energy(N, rodBC)

# Set potential energies
rodBC.potential_energy = (1/2)*kT*phi**2;

# Get (unforced) Lagrangian of the system
L = mech.Lagrangian( N, rodBC );

# Get equation of motion
LM = mech.LagrangesMethod(L, [rodBC], frame=N);
LM.form_lagranges_equations()