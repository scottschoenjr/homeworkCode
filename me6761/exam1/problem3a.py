# -*- coding: utf-8 -*-
"""
Acoustics II Take-Home Exam Problem 3a
"""

import numpy as np
import matplotlib.pyplot as plt

# Create spacial variables
d = 1; # [m]
nPoints = 1000;
x_1d = np.linspace(0, 2*d, 2*nPoints); # [m]
z_1d_above = np.linspace(0, d, nPoints); # [m]
z_1d_below = np.linspace(-d, 0, nPoints); # [m]
x, z_a = np.meshgrid(x_1d, z_1d_above);
x, z_b = np.meshgrid(x_1d, z_1d_below);

# Set material parameters
f = 2000; # Frequency [Hz]
c0 = 100; # Sound speed in fluid [m/s]
rho0 = 1.225; # Density of fluid [kg/m^3]
L = 3E-3; # Plate thickness [m]
U0 = 1E-3; # Oscillation amplitude [m]

# Set panel parameters
rho2 = 2.7E3; # Density of panel [kg/m^3]
E = 70E9; # Elasitc modulus [Pa]
nu = 0.334; # Poisson ratio

# Compute constants
omega = 2*np.pi*f;

# Set up figure
plt.figure(facecolor='white');

# Compute plate velocity 
beta = (E*(L**3))/(12*(1 - nu**2));    
m = rho2/L; # Acoustic mass of panel [kg/m^2]
cpl = ( ((omega**2)*beta)/(rho2*L) )**(0.25);

# Define shorthands
gamma = ( 1 - c0**2/cpl**2  )**0.5;
A = -1j*omega*rho0*c0*U0;
phi_a = 1j*omega*( x/cpl + z_a*gamma/c0 );2
phi_b = 1j*omega*( x/cpl - z_b*gamma/c0 );

# Find pressure fields
p_a = A*np.exp( phi_a );
p_b = A*np.exp( phi_b );

# Normalize
max_a = np.max( np.max( np.abs(p_a) ) );
max_b = np.max( np.max( np.abs(p_b) ) );
normFactor = 1/np.max( [max_a, max_b] );
p_a_Norm = normFactor*p_a;
p_b_Norm = normFactor*p_b;

# Compute velocity fields
vx_a = p_a/( rho0*cpl );
vx_b = p_b/( rho0*cpl );
vz_a = ( gamma/(rho0*c0) )*p_a;
vz_b = -( gamma/(rho0*c0) )*p_b;

# Normalize
max_vxa = np.max( np.max( np.abs(vx_a) ) );
max_vxb = np.max( np.max( np.abs(vx_b) ) );
max_vza = np.max( np.max( np.abs(vz_a) ) );
max_vzb = np.max( np.max( np.abs(vz_b) ) );

normFactor = 1/np.max( [max_vxa, max_vxb, max_vza, max_vzb] );
vx_a_Norm = normFactor*vx_a;
vx_b_Norm = normFactor*vx_b;
vz_a_Norm = normFactor*vz_a;
vz_b_Norm = normFactor*vz_b;

# Plotting options
step = round( nPoints/72 );
scaleFactor = 500;
width = 0.2;

# Plot pressure
colorPlotHandle1 = plt.pcolormesh(x, z_a, np.real(p_a_Norm), \
    cmap='binary',linewidth=0);
colorPlotHandle2 = plt.pcolormesh(x, z_b, np.real(p_b_Norm), \
    cmap='binary',linewidth=0);
colorPlotHandle1.set_edgecolor('face')
colorPlotHandle2.set_edgecolor('face')

plt.ylabel(r'$z$ [m]',fontsize=22, family='serif');
plt.xlabel(r'$x$ [m]',fontsize=22, family='serif');

plt.xticks( fontsize=14, family='serif' );
plt.yticks( fontsize=14, family='serif' );

plt.xlim([0, 2*width]);

plt.clim(-1,1);
cbar = plt.colorbar(colorPlotHandle1, ticks=[-1, -0.5, 0, 0.5, 1]);
cbar.set_label(r'Normalized Pressure', fontsize=16, family='serif');

cbar.ax.yaxis.set_ticklabels( \
    [-1, -0.5, 0, 0.5, 1], \
    family='serif', fontsize=14);

plt.show();

plt.savefig('cbar.png', dpi = 600)

# Plot velocity
plt.figure()

quiverPlotHandle1 = \
    plt.quiver(x[::step,::step], z_a[::step,::step], \
    np.real(vx_a[::step,::step]), np.real(vz_a[::step,::step]), \
    angles='xy', scale=scaleFactor);
quiverPlotHandle2 = \
    plt.quiver(x[::step,::step], z_b[::step,::step], \
    np.real(vx_b[::step,::step]), np.real(vz_b[::step,::step]), \
    angles='xy', scale=scaleFactor);

plt.ylim([-width, width]);
plt.ylabel(r'$z$ [m]',fontsize=22, family='serif');
plt.yticks( fontsize=12, family='serif' );

plt.xlim([0, 2*width]);
plt.xlabel(r'$x$ [m]',fontsize=22, family='serif');
plt.xticks( fontsize=12, family='serif' );

plt.show()


