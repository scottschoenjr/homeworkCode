# -*- coding: utf-8 -*-
"""
Math 6702 Homework 3 Problem 9.6.811
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
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
r = 5;
phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
x = r*np.sin(phi)*np.cos(theta)
y = r*np.sin(phi)*np.sin(theta)
z = r*np.cos(phi)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(facecolor='white')
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface( \
    x, y, z,  rstride=1, cstride=1, color='k', alpha=0.2, linewidth=0.1);
    
# Add in gradient vector
x0 = 3;
y0 = 4;
z0 = 0;
x1 = x0 + 6/5;
y1 = y0 + 8/5;
z1 = z0;
a = Arrow3D([x0, x1], [y0, y1], [z0, z1], mutation_scale=20, 
                lw=3, arrowstyle="-|>", color="k")
ax.add_artist(a)

ax.text( x1, y1, z1 + 1, r'$\nabla F$', None );
    
plt.ylabel(r'$y$',fontsize=26, family='serif',labelpad=20);
plt.yticks( fontsize=16, family='serif' );

plt.xlabel(r'$x$',fontsize=26, family='serif',labelpad=20);
plt.xticks( fontsize=16, family='serif' );

ax.set_zlabel(r'$z$',fontsize=26, family='serif',labelpad=20);
ax.set_zticks( fontsize=16, family='serif' );


