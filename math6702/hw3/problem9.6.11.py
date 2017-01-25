# -*- coding: utf-8 -*-
"""
Math 6702 Homework 3 Problem 9.6.811
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc

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

ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='k', alpha=0.2, linewidth=0.1);
    
plt.ylabel(r'$y$',fontsize=26, family='serif',labelpad=20);
plt.yticks( fontsize=16, family='serif' );

plt.xlabel(r'$x$',fontsize=26, family='serif',labelpad=20);
plt.xticks( fontsize=16, family='serif' );

ax.set_zlabel(r'$z$',fontsize=26, family='serif',labelpad=20);
ax.set_zticks( fontsize=16 );


