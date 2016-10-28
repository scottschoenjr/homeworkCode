# -*- coding: utf-8 -*-
"""
Acoustics II Take-Home Exam Problem 2
"""
import numpy as np
import matplotlib.pyplot as plt

# Define medium properties
rho0 = 1000; # [m/kg^3]
c0 = 1500; # [m/s]

rho1 = 800; # [m/kg^3]
c1 = 1200; # [m/s]
a = 1; # [m]

# Compute parameters
R0 = rho0*c0;
R1 = rho1*c1;
omegaVector = np.logspace(0, 5, 10000);
k0 = omegaVector/c0;
jk0a = 1j*k0*a;
k1 = omegaVector/c1;
jk1a = 1j*k1*a;

# Calculate reflection coefficient
Rnum = 1 - (R0/R1)*( (1 + 1/jk1a)/(1 + 1/jk0a) );
Rden = (R0/R1)*( (1 + 1/jk1a)/(1 + 1/jk0a) ) + ( (1 - 1/jk0a)/(1 + 1/jk0a) );
R = (Rnum/Rden)*np.exp(2*jk0a);
T = (1 + R)*np.exp( 2*jk0a )*(np.exp(-jk0a)/np.exp(-jk1a));

# Plot Reflection and transmission coefficient magnitudes...
plt.figure()
plt.semilogx( k0*a, np.real(R), 'k', label=r'$|\mathcal{R}\,|$' );
plt.semilogx( k0*a, np.abs(T), '--k', label=r'$|\mathcal{T}\,\,|$' );

plt.xlabel(r'$k_{0}a$', fontsize=18, family='serif');
xTickValues = [0.01, 0.1, 1, 10];
xTickLabels = [r'$0.01$', r'$0.1$', r'$1$', r'$10$'];
plt.xticks( xTickValues, xTickLabels, fontsize=13, family='serif' );
plt.xlim([0.01, 20]);

plt.ylabel('Coefficient Magnitude',fontsize=14, family='serif');
plt.yticks( fontsize=12, family='serif' );
plt.ylim([0, 2]);

plt.legend(loc='best', frameon=False, fontsize=12, \
    prop={'family':'serif'});

plt.show()

# ...and Phasses
plt.figure()
plt.semilogx( k0*a, np.angle(R)*180/np.pi, 'k', \
    label=r'$\angle\,\mathcal{R}$' );
plt.semilogx( k0*a, np.angle(T)*180/np.pi, '--k', \
    label=r'$\angle\,\mathcal{T}$' );

plt.xlabel(r'$k_{0}a$', fontsize=18, family='serif');
xTickValues = [0.01, 0.1, 1, 10];
xTickLabels = [r'$0.01$', r'$0.1$', r'$1$', r'$10$'];
plt.xticks( xTickValues, xTickLabels, fontsize=13, family='serif' );
plt.xlim([0.01, 20]);

plt.ylabel('Phase [deg]',fontsize=14, family='serif');
yTickValues = [-180, -90, 0, 90, 180];
plt.yticks( yTickValues, fontsize=12, family='serif' );
plt.ylim([-180, 180]);

plt.legend(loc='best', frameon=False, fontsize=12, \
    prop={'family':'serif'});

plt.show()




