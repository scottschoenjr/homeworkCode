# -*- coding: utf-8 -*-
"""
Problem 3 Homework 11 ME6441
Exercize 7.26 in Ginsberg's Advanced Dynamics
"""

import sympy
import sympy.physics.mechanics as mech

# Set dynamics variables
xc = mech.dynamicsymbols('xc'); # Cart position
xb = mech.dynamicsymbols('xb'); # Box position
xcDot = mech.dynamicsymbols('xcDot'); # Cart velocity
xbDot = mech.dynamicsymbols('xbDot'); # Box velocity

# Create symboloic variables we'll need
m, L, k, g, theta = sympy.symbols('m L k g theta');

# Create reference frame for analysis
N = mech.ReferenceFrame( 'N' );

# Set cart position, mass, and velocity
cPoint = mech.Point( 'cPoint' );
cPoint.set_vel(N, xcDot*N.x);
C = mech.Particle('C', cPoint, 4*m );

# Set box position, mass, and velocity
# Create rotated frame
boxFrame = mech.ReferenceFrame('boxFrame');
boxFrame.orient(N, 'Body', [0, 0, -theta], 'XYZ');
bPoint = mech.Point('bPoint');
bPoint.set_vel(N, xbDot*boxFrame.x);
B = mech.Particle('B', bPoint, m );

# Get the kinetic energy of each component to check results
mech.kinetic_energy(N, C)
mech.kinetic_energy(N, B)

# Set potential energies
C.potential_energy = (1/2)*(2*k)*xc**2;
B.potential_energy = (1/2)*k*xb**2 - m*g*xb*sympy.cos(theta);

# Get (unforced) Lagrangian of the system
L = mech.Lagrangian( N, C, B );

# Get equation of motion
LM = mech.LagrangesMethod(L, [xb, xc], frame=N);
LM.form_lagranges_equations()