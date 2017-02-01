# -*- coding: utf-8 -*-
"""
Math 6702 Homework 3 Problem 9.6.8
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# For drawing vector on 3D plot
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    #

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
    #
#

# Create spatial variables
nPoints = 1000;
x_1d = np.linspace(0.2, 2.9, nPoints);
y_1d = np.linspace(0, 5, nPoints);
x, y = np.meshgrid(x_1d, y_1d);

# Create scalar field
z = (y - 1)/np.sin(x);

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
    
# Add in gradient vector
x0 = np.pi/6;
y0 = 3/2;
z0 = 1;
x1 = x0 + np.sqrt(3);
y1 = y0 + 2;
z1 = z0;
a = Arrow3D([x0, x1], [y0, y1], [z0, z1], mutation_scale=20, 
                lw=3, arrowstyle="-|>", color="k")
ax.add_artist(a)

ax.text( x1 - 1.5, y1-1, z1, r'$\nabla f$', None, fontsize=20 );

# Plot level curve in x-y plane
xLong = np.linspace(0, 3.5, 1000);
linePlotHandle = ax.plot( xLong, np.sin(xLong) + 1, 'k' );

plt.ylim([0, 5]);
plt.ylabel(r'$y$',fontsize=26, family='serif',labelpad=20);
plt.yticks( fontsize=14, family='serif' );

plt.xlim([0, np.pi]);
plt.xlabel(r'$x$',fontsize=26, family='serif',labelpad=20);
plt.xticks( [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], \
    [r'0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', '$\pi$'], \
    fontsize=19, family='serif' );

ax.set_zlabel(r'$z$',fontsize=26, family='serif',labelpad=20);

plt.show()
