# -*- coding: utf-8 -*-
"""
Problem 2 of Homework 4 for Acoustics 2
"""

import numpy as np
import matplotlib.pyplot as plt

# Set material parameters
rho1 = 1000;
c1 = 1500;
Z1 = rho1*c1;

rho2 = 1000;
c2 = 1.1*c1;
Z2 = 1.1*Z1;

# Create incident angle vector
theta_1 = np.linspace(0, np.pi/2, 1000);

# Get transmitted angle
cosThetaT = ( 1 - ((c2/c1)**2)*(np.sin(theta_1))**2 )**(0.5);

# Define reflection and transmission coefficiensts
Rnum = Z2*np.cos(theta_1) - Z1*cosThetaT;
Rden = Z2*np.cos(theta_1) + Z1*cosThetaT;
R = Rnum/Rden;

Tnum = 2*Z2*np.cos(theta_1);
Tden = Z2*np.cos(theta_1) + Z1*cosThetaT;
T = Tnum/Tden;

# Compute intensity transmission coefficient
tau = (Z1/Z2)*(T*np.conjugate(T));
r = R*np.conjugate(R)

# Plot Transmission
plt.figure(facecolor='white');
xTickVals = [0, np.pi/6, np.pi/3, np.pi/2];

plt.subplot(2, 1, 1);
plt.plot(theta_1, abs(T), 'k');

plt.ylim([0, 2]);
plt.ylabel(r'$|\mathcal{T}\,\,|$',fontsize=18, family='serif');
plt.yticks( fontsize=12, family='serif' );

plt.xlim([0, np.pi/2]);
xTickValues = xTickVals;
xTickLabels = [];
plt.xticks( xTickValues, xTickLabels, fontsize=16, family='serif' );

plt.subplot(2, 1, 2);
plt.plot(theta_1, tau, 'k');

plt.ylabel(r'$\tau$',fontsize=18, family='serif');
plt.yticks( fontsize=12, family='serif' );

plt.xlim([0, np.pi/2]);
xTickValues = xTickVals;
xTickLabels = [r'$0$', r'$\pi/6$', r'$\pi/3$', r'$\pi/2$'];
plt.xticks( xTickValues, xTickLabels, fontsize=16, family='serif' );
plt.xlabel(r'Incident Angle $\theta_{1}$ [rad]',fontsize=16, family='serif')

# Plot reflection
plt.figure(facecolor='white');

plt.subplot(2, 1, 1);
plt.plot(theta_1, abs(R), 'k');

plt.ylim([0, 1]);
plt.ylabel(r'$|\mathcal{R}|$',fontsize=18, family='serif');
plt.yticks( fontsize=12, family='serif' );

plt.xlim([0, np.pi/2]);
xTickValues = xTickVals;
xTickLabels = [];
plt.xticks( xTickValues, xTickLabels, fontsize=16, family='serif' );

plt.subplot(2, 1, 2);
plt.plot(theta_1, r, 'k');

plt.ylabel(r'$r$',fontsize=18);
plt.yticks( fontsize=12, family='serif' );

plt.xlim([0, np.pi/2]);
xTickValues = xTickVals;
xTickLabels = [r'$0$', r'$\pi/6$', r'$\pi/3$', r'$\pi/2$']
plt.xticks( xTickValues, xTickLabels, fontsize=16, family='serif' );
plt.xlabel(r'Incident Angle $\theta_{1}$ [rad]',fontsize=16, family='serif');

plt.show()
