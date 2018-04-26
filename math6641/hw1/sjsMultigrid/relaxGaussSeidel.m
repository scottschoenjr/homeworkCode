%**************************************************************************
%
% Perform Gauss-Seidel Decomposition
%
%   Function performs a single iteration of the Gauss-Seidel iteration 
%   process. Solves
%                           Au = f
%   by decomposing A as A = U + L and using forward and back substitution.
%   
%
% Inputs
%   L     - Coefficient Matrix
%   u     - Solution vector
%   f     - Right hand side
%   omega - Weighting term for lower diagonal elements
%
% Outputs
%   u - Next iterated guess for u
%
%**************************************************************************

function u = relaxGaussSeidel(L, u, f, omega)

if omega == 1
	LL = tril(L);
	u = LL\(f - triu(L,1)*u);
else
	% For omega=1, this is equal to tril(L)
	LL = diag(diag(L)) + omega*tril(L,-1);
	u = LL\((1-omega)*diag(L).*u + omega*(f - triu(L,1)*u));
end

end
