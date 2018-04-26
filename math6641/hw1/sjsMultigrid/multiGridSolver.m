% *************************************************************************
% Implementation of Mult-Grid Method
%
%   Solves 2D Poisson equation
%                        u_xx + u_yy = -f
%   with vanishing boundary conditions using the multigrid method.
%
%  Original Code (c) 2012 by Christian B. Mendl (see license)
%  Modified by Scott Schoen Jr
%
%   References:
%	  [1] Trottenberg, Oosterlee, Anton Schuller. "Multigrid",
%         Academic Press (2001)
%	  [2] LeVeque "Finite Difference Methods for Ordinary and Partial
%         Differential Equations" SIAM (2007)
%
% *************************************************************************

clear all
close all
clc

% Set parameters
n = 64; % Number of (square) grid points
maxIterations = 4;
numLevels = 3; % Number of levels to multigrid
nu1 = 1; % Number of pre-smoothing operations
nu2 = 1; % Number of post-smoothing operations
noColorMap = 1;
iterativeSolver = 'G-S'; % 'G-S' or 'Jacobi'/'J'

% Plotting options
plotRHS = 1;


% To plot just the mesh without colormap overlay
whiteMap = ...
    [1, 1, 1; ...
    1, 1, 1 ...
    ];

% Create grid
[x,y] = meshgrid((1:(n-1))/n,(1:(n-1))/n);

% Define forcing term
R0 = 0.25;
f = 0.*x;
f( (x - 0.5).^(2) + (y - 0.5).^(2) < R0.^(2) ) = 1;
fColumnVector = reshape( f, [], 1 );

% Plot f if desired
if plotRHS
    figure();
    set( gca, 'Position', [0.15, 0.15, 0.7, 0.7] );
    surf(x,y,f);
    xlabel( '$x$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
    ylabel( '$y$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
    title( '$f(x, y)$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
    if noColorMap
        colormap( whiteMap );
    end
end

% Compute exact solution
disp( 'Computing with Backslash...' );
tic;
A = poissonStencil2D(n);
u_exact = A \ fColumnVector;
toc;

% Plot exact solution
uPlot = reshape(u_exact, n-1, n-1);
figure();
set( gca, 'Position', [0.15, 0.15, 0.7, 0.7] );
surf(x, y, uPlot );
xlabel( '$x$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
ylabel( '$y$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
zlabel( '$u$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
title( '$-\nabla^{2}u = f$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
if noColorMap
    colormap( whiteMap );
end

% Compute multigrid solution
disp( 'Computing multigrid solution...' );
tic;

% Initialize error
err = cell(2,1);
for gamma = 1:4
    
    % Initialize error vector
    errG = zeros(maxIterations+1,1);
    
    % iterations of multigrid cycle
    
    % Initialize u
    u = zeros( (n-1)^2, 1);
    
    % Compute the error residual
    errG(1) = norm(u - u_exact);
    
    % For each iteration
    for iterationCount = 1 : maxIterations
        u = multigridCycle( ...
            n, gamma, u, @poissonStencil2D, fColumnVector, ...
            nu1, nu2, iterativeSolver );
    end
    
    % store error
    err{gamma} = errG;
end
toc;

% Compute solution with G-S
disp( 'Computing with Iterative Method...' );

% Settings
difference = abs( u - u_exact );
tolerance = max( difference(:) );
maxSolverIterations = 1E6;

% Initialize
uOld = 0.*u_exact;
residual = 1E9;
numIterations = 1;

tic;
while residual > tolerance && numIterations < maxSolverIterations
    
    % Get next iterations and error
    uNew = relaxGaussSeidel( A, uOld, fColumnVector, 1 );
    difference = abs( uNew - u_exact );
    residual = max( difference(:) );
    
    % Update u
    uOld = uNew;
    
    % Increment counter
    numIterations = numIterations + 1;
    
end
toc;

%% Plot results!

% Plot solution from multigrid
uPlot = reshape(u, n-1, n-1);
figure();
set( gca, 'Position', [0.15, 0.15, 0.7, 0.7] );
surf(x, y, uPlot );
xlabel( '$x$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
ylabel( '$y$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
zlabel( '$u$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
title( '$-\nabla^{2}u^{*} = f$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
if noColorMap
    colormap( whiteMap );
end

% Plot errors
gsError = reshape(uNew, n-1, n-1) - reshape(u_exact, n-1, n-1);
mgError = reshape(   u, n-1, n-1) - reshape(u_exact, n-1, n-1);
figure();
set( gcf, 'Position', [100, 100, 1400, 500] );
subplot( 1, 2, 1 );
surf(x, y, gsError );
xlabel( '$x$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
ylabel( '$y$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
zlabel( '$\epsilon$', 'FontSize', 40, 'Interpreter', 'LaTeX' );
title( 'Gauss-Seidel' );
if noColorMap
    colormap( whiteMap );
end
subplot( 1, 2, 2 );
surf(x, y, mgError );
xlabel( '$x$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
ylabel( '$y$', 'FontSize', 28, 'Interpreter', 'LaTeX' );
zlabel( '$\epsilon$', 'FontSize', 40, 'Interpreter', 'LaTeX' );
title( 'Multi-Grid' );
if noColorMap
    colormap( whiteMap );
end

% Plot iteration diagram
figure();
pointVector = 1 : ( 2.*numLevels - 1 );

levelVector = 0.*pointVector + numLevels;
for levelCount = 1 : numLevels
    levelVector( levelCount ) = numLevels - (levelCount - 1);
    levelVector( end - levelCount + 1 ) = numLevels - (levelCount - 1);
    
    % Get string for plot labels
    if levelCount == 1
        levelNoString = '$h$~~';
    else
        multiplier = 2^(levelCount - 1);
        levelNoString = ['$', num2str( multiplier ), 'h$~'];
    end
    yTickLabelVec{levelCount} = levelNoString;
    
end
plot( pointVector, levelVector, '-ok', 'MarkerFaceColor', 'k' );

ylim([0.75, numLevels + 0.25]);
xlim( [0.75, 2.*numLevels - 0.75] );
set( gca, ...
    'XTickLabel', '', ...
    'YTick', [1:numLevels], ...
    'YTickLabel', fliplr( yTickLabelVec ), ...
    'TickLabelInterpreter', 'LaTeX', ...
    'FontSize', 22 ....
    );
box off;

% Plot computation times
figure();
hold all;

nVector = [4, 8, 16, 32, 64, 128, 256];
mTimes  = [0.021554, 0.012373, 0.010811, 0.011598, 0.015459, 0.034194, 0.118016];
mgTimes = [0.084579, 0.126737, 0.278549, 0.713474, 2.419526, 8.392994, 32.911884];
gsTimes = [0.002287, 0.002929, 0.007751, 0.09198, 1.429629, 24.45836, 511.77159];

xVec = log(nVector)./log(2);
matlabLine = plot( xVec, log10(mTimes), ':k', 'LineWidth', 2 );
mgLine = plot( xVec, log10(mgTimes), 'k' );
gsLine = plot( xVec, log10(gsTimes), '--k', 'LineWidth', 2 );

ylabel( '$\log_{10}{t}$', 'Interpreter', 'LaTeX', 'FontSize', 24 );
xlabel( '$\log_{2}{N}$', 'Interpreter', 'LaTeX', 'FontSize', 24 );

legH = legend( [ gsLine, mgLine, matlabLine ], ...
    '~~Iterative', '~~Multi-Grid', '~~MATLAB' );
legH.Location = 'NorthWest';
legH.EdgeColor = 'none';


