# -*- coding: utf-8 -*-
"""
Problem 2 of Homework 4 for Acoustics 2
"""

import numpy as np
import matplotlib.pyplot as plt

# Set material parameters
rho1 = 1000;
c1 = 1000;
Z1 = rho1*c1;

rho2 = 1000;
c2 = 1500;
Z2 = rho2*c2;

# Create incident angle vector
theta_1 = np.linspace(0, np.pi/2, 1000);

# Get transmitted angle
zeta = ( 1 - ((c2/c1)**2)*(np.sin(theta_1))**2 )**(0.5);

# Define reflection and transmission coefficiensts
Rnum = Z1/Z2 - zeta/np.cos(theta_1);
Rden = Z1/Z2 + zeta/np.cos(theta_1);
R = Rnum/Rden;
T = 1 + R;

# Compute intensity transmission coefficient
tau = (Z1/Z2)*(T*np.conj(T));

# Plot!
plt.figure(facecolor='white');

plt.subplot(2, 1, 1);
plt.plot(theta_1, abs(T));
plt.ylabel(r'$\mathcal{T}$',fontsize=18)

plt.xlim([0, np.pi/2]);

plt.subplot(2, 1, 2);
plt.plot(theta_1, tau);
plt.ylabel(r'$\tau$',fontsize=18)

plt.xlim([0, np.pi/2]);
xTickValues = [0, np.pi/6, np.pi/3, np.pi/2]
xTickLabels = [r'$0$', r'$\pi/6$', r'$\pi/3$', r'$\pi/2$']
plt.xticks( xTickValues, xTickLabels, fontsize=16 );

plt.show()
