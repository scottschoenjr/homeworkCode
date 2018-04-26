close all; clear all;
n = 21; % grid points, including internal points
c = 0.1; % constant of tolerance
myforce = @(X,Y) -exp( -(X-0.25).^2 - (Y-0.6).^2 );
 
h     = 1/(n-1);        % grid spacing
tol   = c * h^2;        % tolerance
res   = zeros(n,n);     % storage for residual
u     = res;            % storage for solution
uNew  = u;
k     = 0;                     % loop counter
[X,Y] = meshgrid(0:h:1,0:h:1); % coordinates
f     = myforce(X,Y);
normf = norm( f(:),2); % find norm
 
%Now start the iteration solver, stop when
%relative residual < tolerance
figure;
i     = 2:n-1;
j     = 2:n-1;      %the iterations vectorized
done  = false;
 
while ~done
    k = k+1;
 
    uNew(i,j) = (1/4)*( u(i-1,j) + u(i+1,j) + ...
                    u(i,j-1) + u(i,j+1) - h^2 * f(i,j) );
    res(i,j) = f(i,j) - (1/h^2)*( u(i-1,j) + ...
            u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j) );
 
    % uncomment to see it update as it runs
 
    mesh(X,Y,u);  %hold on;
    title(sprintf('solution at iteration %d',k));
    zlim([0,.07]);
    drawnow;
 
    if norm(res(:),2)/normf < tol
         done = true;
    end;
 
    u = uNew;
end
 
mesh(X,Y,u);
title(sprintf('solution at iteration %d',k))