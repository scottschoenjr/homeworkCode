% Problem 4-10 in _Adv. Eng. Math._, 5e (Zill et al., 2014)

clear all
close all
clc

% Set number of terms
nMax = 500;

% Create x vector
x = linspace( -pi, pi, 5000 );

% Define function
nonzeroIndices = find( x > 0 );
f = 0.*x;
f( nonzeroIndices ) = sin( x(nonzeroIndices ) );

% Compute sum
Sn = 1./pi + 0.5.*sin(x);
for nCount = 2:nMax
    
   n = nCount;
   toAdd = (1./pi).*( ...
       ((-1).^(n) + 1 )./(1 - n.^(2)) ...
       ).*cos( n*x );
   % Add to partial sum
   Sn = Sn + toAdd;
    
end

% Plot function and partial sum
figure()
hold all;

functionPlot = plot( x, f, '--k' );
partialSumPlot = plot( x, Sn, 'k' );
xlabel( '$x$', 'FontSize', 22 );
ylabel( '$f(x)$', 'FontSize', 24 );