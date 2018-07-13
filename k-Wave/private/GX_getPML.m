function D = GX_getPML(Nx, Ny, dx, dt, c, pml_size, pml_alpha)

R = 1e-6;

pml_size = 40;
% pml_alpha = 0.2;
m = 4;

x = 1 : pml_size;

pml_edge1 = (3*c)/(2*pml_size*dx) * ((x-pml_size-1) ./ (0-pml_size)).^m * log10(1/R);
pml_edge2 = (3*c)/(2*pml_size*dx) * (x./pml_size).^m * log10(1/R);

% pml_edge1  = pml_alpha * (c / dx) * ( (x - pml_size - 1) ./ (0 - pml_size) ).^m;
% pml_edge2 = pml_alpha * (c / dx) * ( x ./ pml_size ).^m;


pml_corner = ones(pml_size, pml_size);
for i = 1 : pml_size
    for j = 1 : pml_size
        pml_corner(i, j) = (sqrt((pml_size+1-i).^2 + (pml_size+1-j).^2) / pml_size) .^m;
    end
end
pml_corner = pml_corner * (3*c)/(2*pml_size*dx) * log10(1/R);
% pml_corner = pml_corner * pml_alpha * (c / dx);

D = [pml_corner, pml_edge1' * ones(1, Ny-2*pml_size), fliplr(pml_corner);
    ones(Nx-2*pml_size, 1) * pml_edge1, ones(Nx-2*pml_size, Ny-2*pml_size), ones(Nx-2*pml_size, 1) * pml_edge2;
    flipud(pml_corner), pml_edge2' * ones(1, Ny-2*pml_size), rot90(pml_corner, 2)];