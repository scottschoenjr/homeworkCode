%**************************************************************************
%
% Create 5-point Stencil
%
%   Function to create coefficient matrix for 5 point 2D Laplacian 
%   operator.
%
% Inputs
%   n - Number of interior grid points
%
% Outputs
%   A - Coefficient matrix (n^2 - 1 by n^2 - 1)
%
%**************************************************************************

function A = poisson_stencil2D(n)

% Set step size
h = 1/n;

% example for n = 6:
%	0 1 0 0 0
%	1 0 1 0 0
%	0 1 0 1 0
%	0 0 1 0 1
%	0 0 0 1 0
e = ones(n-1,1);
nD = spdiags([e e],[-1,1],n-1,n-1);

% Kronecker product
nD = kron(speye(n-1),nD) + kron(nD,speye(n-1));

% (n-1)^2 x (n-1)^2 matrix
A = (4*speye((n-1)^2) - nD) / h^2;

end
