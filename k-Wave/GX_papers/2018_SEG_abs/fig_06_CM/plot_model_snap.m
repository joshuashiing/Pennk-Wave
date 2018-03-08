function plot_model_snap(c, p, n_PML, h, scale, clim, figureid, subid)

wf = c + p((n_PML+1) : (end-n_PML), (n_PML+1) : (end-n_PML)) * scale;
[Nx, Ny] = size(wf);
x_axis = (0 : (Nx - 1)) * h;
y_axis = (0 : (Ny - 1)) * h;

figure(figureid)
subplot(subid)
imagesc(y_axis, x_axis, wf, clim);