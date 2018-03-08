function [nx, ny, x, y] = close_grid_2d(kgrid, x0, y0)

% On a 2-D plane, find the grid point closest to (x0, y0)
n = size(x0, 1);
nx = zeros(size(x0));
ny = nx;
x = nx;
y = nx;
for i = 1 : n
    x0i = x0(i);
    y0i = y0(i);
    [xi, nxi] = min(abs(kgrid.x_vec - x0i));
    [yi, nyi] = min(abs(kgrid.y_vec - y0i));
    nx(i) = nxi;
    ny(i) = nyi;
    x(i) = xi;
    y(i) = yi;
end
