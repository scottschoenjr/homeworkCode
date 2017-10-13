# -*- coding: utf-8 -*-
"""
Math 6702 Homework 10 Problem 12.3.47
"""

import numpy as np
import matplotlib.pyplot as plt

# Define forcing function
def squareWave(timeVector, period):
    
    # Initialize
    waveVector= 0*timeVector;
    nPoints = np.size( timeVector )
    
    for timeCount in range(0, nPoints - 1):
        currentTime = timeVector[timeCount];
        
        # Get position within the cycle
        timeWithinCycle = currentTime - np.floor( currentTime/period )*period;
        
        # Determine the fraction of the cycle this is
        fractionOfCycle = timeWithinCycle/period;
        
        if fractionOfCycle < 0.5:
            waveVector[timeCount] = 1;
        else:
            waveVector[timeCount] = -1;
        #
    #
    
    return waveVector
#

# Set time and space dimensions
T = 10; # Square wave period
L = 2*np.pi; # Domain dimension
nValues = [1,2,4,9]

# Create spatial and time variables
numPoints = 1000;
x = np.linspace(0, L, numPoints);
t = np.linspace(0, 3*T, numPoints);

# Compute series solutions for each of the deisred number of terms
numNValues = np.size(nValues);
solutions = np.zeros( (numNValues, numPoints) )
for nValueCount in range(0, numNValues - 1):
    
    # Get number of terms to retain    
    N = nValues[nValueCount];
    
    # Initialize vector to hold current solution
    xn = np.zeros( (1, numPoints) )
    
    # Compute that many coefficients
    for nCount in range(1, N):
        
        # Get current n and coefficients
        n = nCount
        bn = (10/(n*np.pi))*(1 - (-1)**n);
        Bn = bn/( 10 - n**2 );
        
        # Add in contribution
        xn = xn + Bn*(np.sin(n*t) - np.sin(np.sqrt(10)*t)/(n*np.sqrt(10)) );
    #
    
    # Save solution with this number of terms
    solutions[ nValueCount, : ] = xn;

#

# Plot forcing function
plt.figure(facecolor='white');

plotHandle = plt.plot(t, 5*squareWave(t,T), 'k' );

plt.ylim([-5.5, 5.5]);
plt.ylabel(r'$f(t)$',fontsize=22, family='serif');
plt.yticks( fontsize=16, family='serif' );

plt.xlim([0, 2*T]);
plt.xlabel(r'Time [s]',fontsize=16, family='serif');
plt.xticks( fontsize=16, family='serif' );

plt.show()

# Plot solutions with each number of terms
plt.figure(facecolor='white');

for nValueCount in range(0, numNValues - 1):
    # Get number of terms to retain    
    solsPlot = plt.plot( t, solutions[nValueCount, :]);
#

plt.ylim([-5.5, 5.5]);
plt.ylabel(r'$x(t)$',fontsize=22, family='serif');
plt.yticks( fontsize=16, family='serif' );

plt.xlim([0, 2*T]);
plt.xlabel(r'Time [s]',fontsize=16, family='serif');
plt.xticks( fontsize=16, family='serif' );

plt.show()

