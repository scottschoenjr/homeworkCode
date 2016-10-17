# -*- coding: utf-8 -*-
"""
Problem 4 of Homework 3 for Acoustics 2
"""

import numpy as np
import matplotlib.pyplot as plt

# Set material parameters
fVector = np.logspace(3, 4, 1000); # Frequency rnage [Hz]
theta_i_deg = 45; # Incident angle [deg]

# Set panel parameters
in2m = 0.02541;
LVector = in2m*np.array([0.2092, 0.1345, 0.0673, 0.0359]); # Thickness of plate (<< wavelength) [m]
rho2 = 7.8E3; # Density of panel [kg/m^3]

E = 200E9; # Elasitc modulus [Pa]
nu = 0.29; # Poisson ratio


# Set fluid parameters
rho1 = 1.225; # Fluid density [kg/m^3]
c1 = 340; # Sound speed [m/s]

# Compute constants
omega = 2*np.pi*fVector;
theta_i = np.pi*theta_i_deg/180;  # Convert to [rad]

# Set up figure
plt.figure(facecolor='white');
plt.ylabel(r'TL [dB]',fontsize=16, family='serif');
#plt.ylim([10,25]);
#plt.yticks( [10, 15, 20, 25], fontsize=12, family='serif' );

plt.xlabel(r'$f$ [Hz]',fontsize=18, family='serif');
plt.yticks( fontsize=12, family='serif' );

# Plot for each thickness
for L in LVector:
    
    # Get beta
    beta = (E*(L**3))/(12*(1 - nu**2));    
    
    # Acoustic mass of panel [kg/m^2]
    m = rho2/L; 
    cpl = ( ((omega**2)*beta)/(rho2*L) )**(0.25);
    
    # Compute TL for that thickness    
    term1 = ( ((omega*m)/(2*rho1*c1))*np.cos(theta_i) )**2;
    term2 = ( 1 - ((cpl/c1)**4)*(np.sin(theta_i))**4 )**2;
    TL = np.log10( 1 + term1*term2 );
    
    # Label and plot
    L_in = L/in2m;
    labelString = r'$L = %0.3f$ in' % L_in;
    plt.semilogx( fVector, TL, label=labelString, linewidth=2 );
#

plt.legend(frameon=False, loc='upper left')
plt.show()