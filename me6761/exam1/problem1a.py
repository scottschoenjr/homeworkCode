# -*- coding: utf-8 -*-
"""
Acoustics II Take-Home Exam Problem 1a
"""
import numpy as np
import matplotlib.pyplot as plt

# Set ranges
xPos_1d = np.linspace(-1, 1, 400);
yPos_1d = np.linspace(-1, 1, 400);

# Create meshgrid of position
x, y = np.meshgrid( xPos_1d, yPos_1d );

# Create vectors so we can plot wall
xWall = [0, 0];
yWall = [-100, 100];

# Define source and image positions
kd = 0.5;
d = 0.2;
r0 = [d, 0];
ri = [-d, 0];

# Define field contributions
k = kd/d;
r0 = np.sqrt( (x - r0[0])**2 + (y - r0[1])**2 );
p0 = np.exp( -1j*k*r0 )/( 4*np.pi*r0 );

ri = np.sqrt( (x - ri[0])**2 + (y - ri[1])**2 );
pi = -np.exp( -1j*k*ri )/( 4*np.pi*ri );

# Define total field
ptot = p0 + pi;
pNorm = np.abs(ptot)/(np.max(np.max(np.abs(ptot))));
pNorm_dB = 20*np.log10(pNorm);

# Plot result
figHandle = plt.figure(facecolor="white")

colorPlotHandle = plt.pcolormesh(x/d, y/d, pNorm_dB, \
    cmap='viridis',linewidth=0);
colorPlotHandle.set_edgecolor('face')
plt.plot( xWall, yWall, '--w' )

plt.clim(-60,0);
cbar = plt.colorbar();
cbar.set_label(r'Normalized Level [dB]', fontsize=16, family='serif');

plt.xlabel(r'$x/d$',fontsize=20);
plt.xlim([-5, 5]);
#xTickValues = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
#xTickLabels = [r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
plt.xticks( fontsize=14, family='serif' );

plt.ylabel(r'$y/d$',fontsize=20);
plt.ylim([-5, 5]);
#yTickValues = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
#yTickLabels = [r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
plt.yticks( fontsize=14, family='serif' );

titleString = r'Pressure-Release Boundary $kd = %0.1f$' % (kd);
plt.title( titleString, fontsize=18, family='serif' );

plt.show()


