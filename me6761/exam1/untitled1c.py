# -*- coding: utf-8 -*-
"""
Acoustics II Take-Home Exam Problem 1c
"""

import numpy as np
import matplotlib.pyplot as plt

# Set material parameters
f = 5000; # Frequency [Hz]
theta_i_deg = 45; # Incident angle [deg]

# Compute constants
omega = 2*np.pi*f;
theta_i = np.pi*theta_i_deg/180;  # Convert to [rad]

# Set panel parameters
L = 1E-2; # Thickness of plate (<< wavelength) [m]
rho2 = 7.8E3; # Density of panel [kg/m^3]

E = np.logspace( 8, 12, 1000 ); # Elasitc modulus [Pa]
nu = 0.29; # Poisson ratio

# Set fluid parameters
rho1 = 1.225; # Fluid density [kg/m^3]
c1 = 340; # Sound speed [m/s]

# Compute constants
omega = 2*np.pi*f;
theta_i = np.pi*theta_i_deg/180;  # Convert to [rad]

# Set up figure
plt.figure(facecolor='white');

# Compute the transmission loss 
beta = (E*(L**3))/(12*(1 - nu**2));    
    
# Acoustic mass of panel [kg/m^2]
m = rho2/L; 
cpl = ( ((omega**2)*beta)/(rho2*L) )**(0.25);
    
# Compute TL for that thickness    
T = ( 1 + (1j*omega*m/(2*rho1*c1))*(1 - (cpl*np.sin(theta_i)/c1)**4 ) )**-1;
TL = -10*np.log10( abs(T) );

# Label and plot
plt.semilogx( E/1E9, TL, 'k', linewidth=2 );

plt.ylabel(r'TL [dB]',fontsize=16, family='serif');
plt.ylim([50,90]);
plt.yticks( fontsize=14, family='serif' );

plt.xlabel(r'Elastic Modulus [GPa]',fontsize=18, family='serif');
xTickValues = [1, 10, 100, 500]
xTickLabels = [r'$1$', r'$10$', r'$100$', r'$500$'];
plt.xticks( xTickValues, xTickLabels, fontsize=14, family='serif' );

plt.show()