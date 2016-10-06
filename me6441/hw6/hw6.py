# HW 6
# Scott Schoen Jr 20161006

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define constants
L = 10 # Length of arm [m]
g = 9.81 # Gravitational acceleration [m/s^2]
omega = 2*np.sqrt(g/L) # Driving frequency [rad/s]
u0 = 0.1 # Normalized amplitude of the excitation

# Define the ODE system to solve
def simplePendulum( odeVariables, timeVector ):
    
    # Get initial conditions from input
    theta = odeVariables[0]    # Position [rad]
    thetaDot = odeVariables[1] # Angular velocity [rad/s]
    
    # Return ODE array describing motion
    EOM = [ thetaDot, -(g/L)*np.sin(theta) ]
    return EOM
    
#

def forcedPendulum( odeVariables, timeVector ):
    
    # Get initial conditions from input
    theta = odeVariables[0]    # Position [rad]
    thetaDot = odeVariables[1] # Angular velocity [rad/s]
    
    # Return ODE array describing motion
    EOM1 = thetaDot
    EOM2 = -(g/L)*np.sin(theta) + ( omega**2*u0*np.sin(omega*timeVector) )*np.sin(theta)
    EOM = [ EOM1, EOM2  ]
    return EOM
    
#
    
# Set initial conditions
theta0 = 1; # [rad]
thetaDot0 = 0.5; # [rad/s] 
initialConditions = [ theta0, thetaDot0 ]

# Set time limits
tMin = 0  # [s] Time of theta0 and thetaDot0
tMax = 20 # [s]
tVector = np.linspace( tMin, tMax, 200 )

# Get solution
simpleSolution = odeint( simplePendulum, initialConditions, tVector )
solution = odeint( forcedPendulum, initialConditions, tVector )
theta = solution[:,0]
thetaDot = solution[:,1]
    
# Plot results
plt.figure(facecolor='white')
plt.plot( tVector, theta, 'k' )
plt.xlabel( 'Time [s]', fontsize=16 )
plt.ylabel( r"Position $\theta$ [rad]", fontsize=16 )
plt.show()

plt.figure(facecolor='white')
plt.plot( tVector, thetaDot, 'k' )
plt.xlabel( 'Time [s]', fontsize=16 )
plt.ylabel( r"Angular Velocity $\dot{\theta}$ [rad/s]", fontsize=16 )
plt.show()