# HW 6
# Scott Schoen Jr 20161006

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
from scipy.integrate import odeint

# Define constants
L = 1 # Length of arm [m]
g = 9.81 # Gravitational acceleration [m/s^2]

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
    EOM2 = -(g/L)*np.sin(theta) \
            + ( omega**2*u0*np.sin(omega*timeVector) )*np.sin(theta)
    EOM = [ EOM1, EOM2  ]
    return EOM
    
#
    
# ------- Problem 1(c) ---------
u0 = 0.1; # Amplitude relative to length of arm L
omega = 2*np.sqrt( g/L );
    
# Set initial conditions
theta0 = 0.5; # [rad]
thetaDot0 = 0.0; # [rad/s] 
initialConditions = [ theta0, thetaDot0 ]

# Set reference initial conditions
theta0Ref = 0.5; # [rad]
thetaDot0Ref = 0.0; # [rad/s] 
initialConditionsRef = [ theta0Ref, thetaDot0Ref ]

# Set time limits
tMin = 0   # [s] Time of theta0 and thetaDot0
tMax = 100 # [s]
dt = 0.01; # Time Step [s]
nPoints = np.round( (tMax - tMin)/dt );
tVector = np.linspace( tMin, tMax, nPoints )

# Get solutions
solution = odeint( forcedPendulum, \
      initialConditions, tVector )
theta = solution[:,0]
thetaDot = solution[:,1]
    
refSolution = odeint( simplePendulum, \
      initialConditionsRef, tVector )
thetaRef = refSolution[:,0]
thetaDotRef = refSolution[:,1]

# Plot results
fig = plt.figure(facecolor='white');
ax = fig.add_subplot(111);
solutionPlot, = plt.plot( tVector, theta, 'k', \
     label=r'$u_{0} = 0.1$' );
referencePlot, = plt.plot( tVector, thetaRef, '--b', \
     label=r'$u_{0} = 0$', linewidth=1 );
     
# Set x-range of interest
x0 = 0;
x1 = 20; 

# Formatting...
plt.xlabel( 'Time [s]', fontsize=16, family='serif' )
plt.ylabel( r"Position $\theta$ [rad]", fontsize=16, family='serif' )
rc('font', family='serif');
yTickValues = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
yTickLabels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
plt.yticks( yTickValues, yTickLabels, fontsize=16 );
# Plot shaded region
y0 = min(yTickValues);
width = x1 - x0;
height = max(yTickValues) - y0;
rect1 = mpl.patches.Rectangle((x0,y0), width, height, \
     color='#000000', alpha=0.2);
ax.add_patch(rect1)
plt.legend(frameon=False)
plt.show()

# Plot inset
plt.figure(facecolor='white')
solutionPlot, = plt.plot( tVector, theta, 'k', \
     label=r'$u_{0} = 0.5$ rad/s' );
referencePlot, = plt.plot( tVector, thetaRef, '--b', \
     label=r'$u_{0} = 0$ rad/s', linewidth=1 );

# Formatting...
plt.xlabel( 'Time [s]', fontsize=16, family='serif' )
plt.ylabel( r"Position $\theta$ [rad]", fontsize=16, family='serif' )
rc('font', family='serif');
plt.xlim([x0, x1]);
yTickValues = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
yTickLabels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
plt.yticks( yTickValues, yTickLabels, fontsize=16 );
plt.show()

# ------- Problem 1(d) ---------
u0 = 0.1; # Amplitude relative to length of arm L
omega = 40*np.sqrt( g/L );
    
# Set initial conditions
theta0 = 0.96*np.pi; # [rad]
thetaDot0 = 0.5; # [rad/s] 
initialConditions = [ theta0, thetaDot0 ]

# Set reference initial conditions
theta0Ref = 0.5; # [rad]
thetaDot0Ref = 0.0; # [rad/s] 
initialConditionsRef = [ theta0Ref, thetaDot0Ref ]

# Set time limits
tMin = 0   # [s] Time of theta0 and thetaDot0
tMax = 10 # [s]
dt = 0.01; # Time Step [s]
nPoints = np.round( (tMax - tMin)/dt );
tVector = np.linspace( tMin, tMax, nPoints )

# Get solution
solution = odeint( forcedPendulum, \
      initialConditions, tVector )
theta = solution[:,0]
thetaDot = solution[:,1]


# Plot results
plt.figure(facecolor='white')
solutionPlot, = plt.plot( tVector, theta, 'k' );

# Formatting...
plt.xlabel( 'Time [s]', fontsize=16, family='serif' )
plt.ylabel( r"Position $\theta$ [rad]", fontsize=16, family='serif' )
rc('font', family='serif');
yTickValues = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
yTickLabels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$',]
plt.yticks( yTickValues, yTickLabels, fontsize=16 );
plt.show()

