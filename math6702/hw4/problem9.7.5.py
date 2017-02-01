# -*- coding: utf-8 -*-
"""
Math 6702 Homework 3 Problem 9.7.5
"""

import numpy as np
import matplotlib.pyplot as plt

# Create spatial variables
nPoints = 1000;
x_1d = np.linspace(-1, 1, nPoints);
y_1d = np.linspace(-1, 1, nPoints);
x, y = np.meshgrid(x_1d, y_1d);

# Create vector field
Fx = 0*x;
Fy = 1*y;

# Set up figure
plt.figure(facecolor='white');

step = round( nPoints/14 );
scaleFactor = 15;
width = 0.2;

quiverPlotHandle = \
    plt.quiver(x[::step,::step], y[::step,::step], \
    Fx[::step,::step], Fy[::step,::step], \
    angles='xy', scale=scaleFactor);

plt.ylim([-1, 1]);
plt.ylabel(r'$y$',fontsize=22, family='serif');
plt.yticks( fontsize=16, family='serif' );

plt.xlim([-1, 1]);
plt.xlabel(r'$x$',fontsize=22, family='serif');
plt.xticks( fontsize=16, family='serif' );

plt.show()

