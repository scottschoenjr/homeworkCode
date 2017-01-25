# -*- coding: utf-8 -*-
"""
Math 6702 Homework 3 Problem 9.6.8
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

# Create spatial variables
nPoints = 1000;
x_1d = np.linspace(-5, 5, nPoints);
y_1d = np.linspace(-3, 3, nPoints);
x, y = np.meshgrid(x_1d, y_1d);

# Create scalar field
z = (y - 1)/np.sin(x);
# z[ z > 2 ] = np.NAN;
# z[ z < 0 ] = np.NAN;

# Compute level curve value
x0 = np.pi/6;
y0 = 3/2;
c = (y0 - 1)/np.sin(x0);

# Set up figure
fig = plt.figure(facecolor='white');
ax = fig.gca(projection='3d')

tol = 0.0005;
levels = [ c - tol, c + tol ];

maxValue = 2;

surfacePlotHandle = ax.plot_surface(x, y, z, \
    cmap='viridis', linewidth=0, antialiased=False );

    
# contourPlotHandle = \
#    plt.contour(x, y, z, levels, colors='k');
    
    
ax.set_zlim(0, 2);

# fig.colorbar(surfacePlotHandle);
v = np.linspace( 0, 2, 32);
x = plt.colorbar(surfacePlotHandle, ticks=v)

plt.ylim([-0.05, 2.05]);
plt.ylabel(r'$y$',fontsize=22, family='serif');
plt.yticks( fontsize=16, family='serif' );

plt.xlim([0, np.pi]);
plt.xlabel(r'$x$',fontsize=22, family='serif');
plt.xticks( fontsize=16, family='serif' );

plt.show()
