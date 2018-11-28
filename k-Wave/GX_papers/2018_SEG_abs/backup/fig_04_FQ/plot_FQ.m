function plot_FQ(case_list, h, figureid)

casename1 = case_list{1};
load(casename1);
p1 = p;
[Nx, Ny] = size(p1);
p_FQ = zeros(size(p1));
NQx = floor(Nx / 2);
NQy = floor(Ny / 2);
p_FQ(1:NQx, 1:NQy) = p1(1:NQx, 1:NQy);

casename2 = case_list{2};
load(casename2);
p2 = p;
p_FQ(1:NQx, (NQy+1):Ny) = p2(1:NQx, (NQy+1):Ny);

casename3 = case_list{3};
load(casename3);
p3 = p;
p_FQ((NQx+1):Nx, 1:NQy) = p3((NQx+1):Nx, 1:NQy);

casename4 = case_list{4};
load(casename4);
p4 = p;
p_FQ((NQx+1):Nx, (NQy+1):Ny) = p4((NQx+1):Nx, (NQy+1):Ny);

x_axis = (0 : (Nx - 1)) * h;
y_axis = (0 : (Ny - 1)) * h;
u_lim = max(abs(p_FQ(:)));
figure(figureid);
imagesc(y_axis, x_axis, p_FQ, [-u_lim, u_lim]);
colorbar;