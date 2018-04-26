%**************************************************************************
%
% Get Nonlinear Matrices
%
%   Function computes the nonlinear "C" matrices for the FEM problem.
%
%**************************************************************************

function [ C1, C2, C3, gC1, gC2, gC3 ] = ...
    getNonlinearMatrices( Mesh, n, rho0, c0, beta )

% Function that build the nonlinear matrices using element matrices 
% from Eq. (4.19) [Dirske 2014]
C1EI = cell(n);
C2EI = cell(n);
C3EI = cell(n);

% Nonlinear factor shorthand
nlFactor = 2.*beta./( rho0.*c0.^(4) );

for i=1:n
    C1EI{i} = nlFactor * 1/12 * [0 1; 0 1] * (Mesh(i+1) - Mesh(i));
    C2EI{i} = nlFactor * 1/12 * [3 1; 1 3] * (Mesh(i+1) - Mesh(i));
    C3EI{i} = nlFactor * 1/12 * [1 0; 1 0] * (Mesh(i+1) - Mesh(i));
end

C1EI{1}(1,2) = 0;
C3EI{n}(2,1) = 0;

% Create sparse matrices
C1 = spalloc( n+1, n+1, 3*(n+1) );
C2 = spalloc( n+1, n+1, 3*(n+1) );
C3 = spalloc( n+1, n+1, 3*(n+1) );

% Fill in elements
for i=1:n
    
    C1(i,i) = C1(i,i) + C1EI{i}(1,1);
    C1(i,i+1) = C1(i,i+1) + C1EI{i}(1,2);
    C1(i+1,i) = C1(i+1,i) + C1EI{i}(2,1);
    C1(i+1,i+1) = C1(i+1,i+1) + C1EI{i}(2,2);
    
    C2(i,i) = C2(i,i) + C2EI{i}(1,1);
    C2(i,i+1) = C2(i,i+1) + C2EI{i}(1,2);
    C2(i+1,i) = C2(i+1,i) + C2EI{i}(2,1);
    C2(i+1,i+1) = C2(i+1,i+1) + C2EI{i}(2,2);
    
    C3(i,i) = C3(i,i) + C3EI{i}(1,1);
    C3(i,i+1) = C3(i,i+1) + C3EI{i}(1,2);
    C3(i+1,i) = C3(i+1,i) + C3EI{i}(2,1);
    C3(i+1,i+1) = C3(i+1,i+1) + C3EI{i}(2,2);
    
end

% Define RHS vector
gC1 = C1(2,1);
gC2 = C2(2,1);
gC3 = C3(2,1);

C1 = C1(2:n,2:n);
C2 = C2(2:n,2:n);
C3 = C3(2:n,2:n);

end
