%**************************************************************************
%
% Perform Multi-Grid Cycle
%
%  Following procedure from p. 107 of Ref. 1
%
% Inputs
%   n              - Number of interior grid points
%   gamma          - Maximum number of grid reductions (coarsest grid will
%                    have N/2^gamma interior points
%   u              - Solution vector (initial guess)
%   Lop            - Function handle for Laplacian stencil
%   f              - RHS of equation
%   preIterations  - Number of initial smoothing operations
%   postIterations - Number of final smoothing operations         
%   smoothing      - Smoothing technique
%                    'Gauss-Seidel'/'G-S'/'G' or 'Jacobi'/'J'
%
% Outputs
%   u - Corrected solution from multi-grid cycle
%
%**************************************************************************

function u = multigridCycle(n, gamma, u, Lop, f, ...
    preIterations, postIterations, smoothing)

% Make sure n is a powre of 2
if mod(n,2)~=0
	error('n must be a power of 2.');
end

% Get coefficient matrix
A = Lop(n);

% Set smoothing method
switch lower( smoothing )
	case {'gaussseidel', 'sauss-seidel', 'g', 'gs', 'g-s'}
		% set omega = 1
		smooth = @(L,u,f) relaxGaussSeidel(L, u, f, 1);
	case {'jacobi', 'j'}
		% for smoothing of Poisson equation, omega = 4/5 is optimal
		smooth = @(L, u, f) relaxJacobi(L,u,f,4/5);
	otherwise
		error('Unknown smoothing method. Specify Jacobi or Gauss-Seidel.');
end

% STEP 1: Apply pre-smoothing
for count = 1:preIterations
	u = smooth(A,u,f);
end

% STEP 2: Compute error residual
rV = f - A*u;

% STEP 3: Coarsen grid
rV = restrictionFW2D(n)*rV;

% STEP 4: compute an (approximate) solution v^tilde_{n-1}
if n==2
	v = Lop(n/2)\rV;
else
	v = zeros((n/2-1)^2,1);	% start with zero
	for j=1:gamma
		v = multigridCycle( ...
              n/2, gamma, v, Lop, rV, ...
              preIterations, postIterations, smoothing ...
            );
	end
end

% STEP 5: Interpolate the error to the (original) finer mesh and apply
% correction
v = interpolation2D(n/2)*v;
u = u + v;

% STEP 6: Smooth out errors in corrected solution
for count = 1:postIterations
	u = smooth(A,u,f);
end

end

% References
%  [1] LeVeque "Finite Difference Methods for Ordinary and Partial
%      Differential Equations" SIAM (2007)