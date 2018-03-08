function f_vec = make_dual_dim(Nt, dt)
% Create a dual dimension for Fourier transform

if rem(Nt, 2) == 0
    f_vec = (-Nt/2 : (Nt/2 - 1))' / Nt;
else
    f_vec = (-(Nt-1)/2 : (Nt-1)/2)' / Nt;
end

f_vec(floor(Nt / 2) + 1) = 0;
f_vec = f_vec ./ dt;