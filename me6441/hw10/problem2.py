# -*- coding: utf-8 -*-
"""
Problem 2 Homework 10 ME6441
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
from scipy.integrate import odeint

# Define parameters
L = 2; # Box length [m]
R = 1; # Semicircle radius [m]
h = 0.04; # Box height [m]
g = 9.81; # Gravitational acceleration [m/s^2]

# Define the ODE system to solve
def rockingBox( odeVariables, timeVector ):
    
    # Get initial conditions from input
    phi = odeVariables[0]    # Position [rad]
    phiDot = odeVariables[1] # Angular velocity [rad/s]
    
    # Return ODE array describing motion
    phiDdotNum = -g*( R*phi*np.cos(phi) - (h/2)*np.sin(phi) ) \
        - ((R**2)*(phiDot**2)*phi);
    phiDdotDen = ( (L**2)/12 + (h**2)/3 + (R**2)*(phi**2) );
    phiDdot = phiDdotNum/phiDdotDen;
    EOM = [ phiDot, phiDdot ];
    return EOM
    
#

# Set initial conditions
tMin = 0; # Start time [s]
tMax = 5;  # End time [s]
phi0 = 1; # Initial angle [rad]
phiDot0 = 0; # Initial angular velocity [rad/s]
initialConditions = [phi0, phiDot0];

# Set time limits
dt = 0.01; # Time Step [s]
nPoints = np.round( (tMax - tMin)/dt );
tVector = np.linspace( tMin, tMax, nPoints )

# Get solutions
solution = odeint( rockingBox, \
      initialConditions, tVector )
phi = solution[:,0]
phiDot = solution[:,1]
phi_deg = 180*phi/np.pi;

# Plot results
fig = plt.figure(facecolor='white');
ax = fig.add_subplot(111);
solutionPlot, = plt.plot( tVector, phi_deg, 'k', \
     label=r'$u_{0} = 0.1$' );
     
# Set range of interest
plt.xlim([0, 5]);
plt.ylim([-90, 90]);

# Formatting...
plt.xlabel( 'Time [s]', fontsize=16, family='serif' )
plt.ylabel( r'$\phi$ [deg]', fontsize=18, family='serif' )
rc('font', family='serif');
yTickValues = [-90, -45, 0, 45, 90]
plt.yticks( yTickValues, fontsize=16 );
