% Create arbitrary vector field

x = linspace( 0.8, 1.5, 6 );
y = linspace( 0.8, 1.5, 6 );
[x, y] = meshgrid( x, y );
% xShift = linspace( 0, 0.7, 6 )';
% x = x + repmat( xShift, 1, 6 )';

u = cos(x);
v = sin(x);

qPlot = quiver( x, y, u, v, 'k' );
qPlot.LineWidth = 2;
qPlot.Color = [0, 0, 0];
qPlot.MaxHeadSize = 0.1;
qPlot.AlignVertexCenters = 'on';
axis off;