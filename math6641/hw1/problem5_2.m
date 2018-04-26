% Problem 5

clear all
close all
clc

% Create grid points
h0 = 1E-2;
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

plotProgress = 1;
maxIterations = 100;

deltaX = h0;
deltaY = h0;

xVec = xMin - deltaX : deltaX : xMax + deltaX;
yVec = yMin - deltaY : deltaY : yMax + deltaY;

[x, y] = meshgrid( xVec, yVec );

% Define forcing and boundary condition functions
f = -1.*exp( -(20./(xMax - xMin)).*(x-1./2).^2 - (20./(yMax - yMin)).*(y-1./3).^2 );

% Initialize U
U = 0.*f;

% Let's keep it simple and say U = 0 on the boundaries
U( :, 1 ) = 0;
U( :, end ) = 0;
U( 1, : ) = 0;
U( end, : ) = 0;

% Create the FD matrix
Nx = length( xVec );
Ny = length( yVec );
j = 2 : Nx - 1;
k = 2 : Ny - 1;

% Set up matrix equation
( U( j + 1, k ) + U( j - 1, k ) ).*deltaY.^(2) + ...
    U(j, k + 1) + U( j, k - 1 )

% Iterate U at current resolution
numIterations = 0;
while numIterations < maxIterations
    
    % Store old U
    Uold = U;
    
    % Update with Jacobi method
    U(j, k) = ...
        ( ...
        ( Uold( j - 1, k ) + Uold( j + 1, k ) ).*deltaX.^(2) + ...
        ( Uold( j, k - 1 ) + Uold( j, k + 1 ) ).*deltaY.^(2) - ...
        deltaX.^(2).*deltaY.^(2).*f( j, k ) ...
        ) ...
        ./ ...
        ( 2.*( deltaX.^(2) + deltaY.^(2) ) );
    
    
    % Plot progress
    if plotProgress
        figure(1)
        mesh( x, y, U );
        shading flat;
        xlabel( 'x' );
        ylabel( 'y' );
        title( ['Iteration ', num2str(numIterations) ] );
                drawnow();
    end
    
    % Increment
    numIterations = numIterations + 1;
    
end

% Plot final solution
mesh(x, y, U);
xlabel( 'x' );
ylabel( 'y' );
zlabel( 'U' );



