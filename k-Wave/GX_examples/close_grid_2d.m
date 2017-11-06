function [nx, ny, x, y] = close_grid_2d(kgrid, x0, y0)

% On a 2-D plane, find the grid point closest to (x0, y0)
[x, nx] = min(abs(kgrid.x_vec - x0));
[y, ny] = min(abs(kgrid.y_vec - y0));