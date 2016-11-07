# -*- coding: utf-8 -*-
"""
ME 6761 Acoustics II, Homework 5, Problem 2c
"""

import numpy as np
import matplotlib.pyplot as plt

# Set distance variable
r = np.logspace( -2, 2.5, 1000 )

# Define source parameters
f = 1E3; # [Hz]
xs = 0.3; # [m]

omega = 2*np.pi*f;
c = 343; # Assume air [m/s]
k = omega/c;

# Define Fraunhoeffer and Fresnel Parameters
fresnelParameter = xs/r;
fraunhoeferParametrer = (k*xs**2)/r;

# Plot each
plt.loglog( r, fresnelParameter, 'k', label=r"$r'/r$" );
plt.loglog( r, fraunhoeferParametrer, '--k', label=r"$kr'\,^{2}/r$" );

plt.xlabel(r'$r$ [m]', fontsize=18, family='serif');
xTickValues = [0.1, 1, 10, 100];
xTickLabels = [r'$0.1$', r'$1$', r'$10$', r'$100$'];
plt.xticks( xTickValues, xTickLabels, fontsize=15, family='serif' );

plt.ylabel('Parameter Value',fontsize=14, family='serif');
yTickValues = [0.1, 1, 10];
yTickLabels = [r'$0.1$', r'$1$', r'$10$'];
plt.yticks( yTickValues, yTickLabels, fontsize=15, family='serif' );

plt.ylim([0.01, 10]);
plt.xlim([0.1, 100]);

plt.legend(loc='best', frameon=False, fontsize=18, \
    prop={'family':'serif'});