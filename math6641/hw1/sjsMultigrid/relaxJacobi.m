%**************************************************************************
%
% Perform Jacobi Iteration
%
%   Function performs a single iteration of the Jacobi iteration 
%   process. Solves
%                           Au = f
%   by decomposing A as A = D + M and using forward and back substitution.
%   
%
% Inputs
%   L     - Coefficient Matrix
%   u     - Solution vector
%   f     - Right hand side
%   omega - Weighting term 
%
% Outputs
%   uNew - Next iterated guess for u
%
%**************************************************************************

function uNew = relaxJacobi(L, u, f, omega)

uNew = u + omega./diag(L).*(f - L*u);

end

% References
%   [1] Ulrich Trottenberg, Cornelius W. Oosterlee, Anton Schuller.
%		Multigrid, Academic Press (2001)
