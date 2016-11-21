# -*- coding: utf-8 -*-
"""
Problem 5 Homework 11 ME6441
"""

import sympy
import sympy.physics.mechanics as mech

# Set dynamics variables
theta = mech.dynamicsymbols('theta');
thetaDot = mech.dynamicsymbols('thetaDot');

# Create symboloic variables we'll need
m, L, k, g = sympy.symbols('m L k g');

# Set velocity of point B
vB = L*thetaDot*sympy.cos(theta)/2; # (y-direction)

# Create reference frame for analysis
N = mech.ReferenceFrame( 'N' );

# Set origin
O = mech.Point('O');
O.set_vel(N, 0*N.x + 0*N.y + 0*N.z);

# Set point B
B = mech.Point( 'B' );
B.set_vel(N, vB*N.y);

# Create vertical rod (this one isn't moving)
I1 = mech.inertia(N, 0, 0, 0);
rod1 = mech.RigidBody( 'rod1', B, N, m, (I1, B) );

# Create rod AB
Izz = m*L**2/12;
I2 = mech.inertia(N, 0, 0, Izz);
rod2Frame = mech.ReferenceFrame('rod2Frame');
rod2Frame.set_ang_vel( N, -thetaDot*N.z );
rod2 = mech.RigidBody( 'rod2', B, rod2Frame, m, (I2, B) );

# Get the kinetic energy of each component tp check results
mech.kinetic_energy(N, rod1)
mech.kinetic_energy(N, rod2)

# Set potential energies
deltaz = (L/2)*( sympy.sin(theta) - 1/2 );
V1 = m*g*deltaz; # Rod 1
V2 = m*g*deltaz; # Rod 2
Vs = (1/2)*k*deltaz**2; # Spring
rod1.potential_energy = V1 + Vs;
rod2.potential_energy = V2;

# Get (unforced) Lagrangian of the system
L = mech.Lagrangian( N, rod1, rod2 );

# Get equation of motion
LM = mech.LagrangesMethod(L, [theta], frame=N);
LM.form_lagranges_equations()
