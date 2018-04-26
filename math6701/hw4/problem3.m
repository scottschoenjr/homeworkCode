% Problem 3

clear all
close all
clc

% Find nth roots of z
z = -2;
n = 4;

r = abs( z );
theta = angle( z );

k = 0 : n - 1;
wk = r.^(1./n).*( ...
    cos( (theta + 2.*pi.*k )./n ) + 1i.*sin( (theta + 2.*pi.*k )./n ) ...
    );

% Plot real and imaginary parts
figure();
hold all;
% axis off;

axis equal;
xlim( r.*[-1.2, 1.3] );
ylim( r.*[-1.2, 1.3] );

% Draw axes
axHandle = gca;
axHandle.XAxisLocation = 'origin';
axHandle.YAxisLocation = 'origin';
plot( max(xlim), 0, 'k>', 'MarkerSize', 6, 'MarkerFaceColor', 'k' );
plot( 0, max(ylim), 'k^', 'MarkerSize', 6, 'MarkerFaceColor', 'k' );

% Label axes
set( gca, ...
    'XTick', [-1, -0.5, 0, 0.5, 1], ...
    'YTick', [-1, -0.5, 0.5, 1], ...
    'YTickLabels', {'$-i\,\,$';'';'';'$i\,\,$'} ...
    );
axHandle.YAxis.FontSize = 24;


% Plot roots
plot( real(wk), imag(wk), 'ko', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k' );
