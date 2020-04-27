% Problem 2

clear all
close all
clc

rVec = [0.5, 1, 1.5, 2, 2.5];
kx = linspace( 0, 4.*pi, 1000 );

figure();
hold all;
set( gca, 'FontSize', 20 );

% Plot magnitude of rho for each r
for rCount = 1 : length(rVec)
   
    r = rVec( rCount );
    rho = cos( kx ) + 1i.*r.*sin( kx );
    
    % Plot
    plot( kx, sqrt( rho.*conj(rho) ), 'k' );
    
    % Create text box to move by hand
    if rCount == 1
    annotation( 'textbox', ...
        'String', ['$r = ', num2str(r), '$'], ...
        'EdgeColor', 'none', ...
        'Interpreter', 'LaTeX', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 17 ...
        );
    else
        annotation( 'textbox', ...
        'String', ['$', num2str(r), '$'], ...
        'EdgeColor', 'none', ...
        'Interpreter', 'LaTeX', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 17 ...
        );
    end
end

xlabel( '$\xi \Delta x$', 'interpreter', 'LaTeX', 'FontSize', 28 );
ylabel( '$|\rho|$', 'interpreter', 'LaTeX', 'FontSize', 32 );

xlim([0, 3.*pi]);
ylim([0, 2.8]);

set( gca, ...
    'Position', [0.15, 0.15, 0.8, 0.8], ...
    'YTick', [0:0.5:2.5], ...
    'XTick', [0, pi, 2.*pi, 3.*pi], ...
    'XTickLabel', {'$0$';'$\pi$';'$2\pi$';'$3\pi$'} ...
    );