# -*- coding: utf-8 -*-
"""
Problem 4 of Homework 3 for Acoustics 2
"""

import numpy as np
import matplotlib.pyplot as plt

# Set material parameters
fVector = 5000; # Frequency [Hz]
theta_i_deg = 45; # Incident angle [deg]

# Set panel parameters
in2m = 0.02541;
L = in2m*np.linspace(0.1, 0.3, 100000); # Thickness of plate (<< wavelength) [m]
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
plt.ylim([40,80]);
#plt.yticks( [10, 15, 20, 25], fontsize=12, family='serif' );

plt.xlabel(r'$L$ [in]',fontsize=18, family='serif');
plt.yticks( fontsize=12, family='serif' );

# Compute the transmission loss 
beta = (E*(L**3))/(12*(1 - nu**2));    
    
# Acoustic mass of panel [kg/m^2]
m = rho2/L; 
cpl = ( ((omega**2)*beta)/(rho2*L) )**(0.25);
    
# Compute TL for that thickness    
T = ( 1 + (1j*omega*m/(2*rho1*c1))*(1 - (cpl*np.sin(theta_i)/c1)**4 ) )**-1;
TL = -10*np.log10( abs(T) );

# Label and plot
L_in = L/in2m;
plt.plot( L_in, TL, 'k', linewidth=2 );

plt.show()