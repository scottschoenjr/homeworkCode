# -*- coding: utf-8 -*-
"""
Problem 4 of Homework 4 for Acoustics 2
"""

import numpy as np
import matplotlib.pyplot as plt

# Set material parameters
fVector = np.logspace(2, 5, 100000); # Frequency rnage [Hz]
incidentAngles = np.array([15, 30, 45]); # Incident angles [deg]

# Set panel parameters
L = 1E-3; # Thickness of plate (<< wavelength) [m]
rho2 = 1.1E3; # Density of panel [kg/m^3]
m = rho2/L; # Acoustic mass of panel [kg/m^2]
T = 10; # Tension

# Set fluid parameters
rho1 = 1.24; # Fluid density [kg/m^3]
c1 = 343; # Sound speed [m/s]

# Compute constants
omega = 2*np.pi*fVector;
cpl = ( ((omega**2)*T)/(rho2*L) )**(0.25);
incidentAngles = np.pi*incidentAngles/180; # Convert to [rad]

# Set up figure
plt.figure(facecolor='white');
plt.ylabel(r'TL [dB]',fontsize=16, family='serif');
plt.ylim([10,25]);
plt.yticks( [10, 15, 20, 25], fontsize=12, family='serif' );

plt.xlabel(r'$f$ [Hz]',fontsize=18, family='serif');
plt.yticks( fontsize=12, family='serif' );


# Compute transmission loss for each incident angle
for theta_i in incidentAngles:
    
    theta_i_deg = round(180*theta_i/np.pi);
    term1 = ( ((omega*m)/(2*rho1*c1))*np.cos(theta_i) )**2;
    term2 = ( 1 - ((cpl/c1)**4)*(np.sin(theta_i))**4 )**2;
    TL = np.log10( 1 + term1*term2 );
    labelString = r'$\theta_{i} = %d ^{\rm o}$' % theta_i_deg;
    plt.semilogx( fVector, TL, label=labelString, linewidth=2 );
#

#plt.legend(frameon=False, loc='upper left')
plt.show()
