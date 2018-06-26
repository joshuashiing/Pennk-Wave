clear;
dx = 0.1;
x = 0 : dx : 100;
x2 = x(1 : end-1) + dx/2;
f = 0.3;
w = 2 * pi * f;
c = 3;
k = w / c;

ki = 1i * 0.1;
k = k + ki;
t = 1;
u = exp(1i * (k * x) - w * t);
% u = real(u);

[k_vec, uu_vec] = comp_spec(x, u);
k_vec = k_vec * 2 * pi;

% [~, ik] = min(abs(k_vec - real(k)));
% k_vec(ik) = k_vec(ik) + ki;

luu_vec = transpose(uu_vec) .* k_vec.^2;
lu = -ifft(ifftshift(luu_vec));


lu_ref = diff(u, 2) / (dx^2);
du_ref = diff(u) / dx;
[k_vec2, luu_vec_ref] = comp_spec(x(2:end-1), lu_ref);
k_vec2 = k_vec2 * 2 * pi;


% plot(x(2:end-1), real(lu_ref), 'k'); hold on;
% plot(x, real(lu), 'r--');
% plot(x, real(u), 'b');
% xlim([10 90]);

plot(x2, real(du_ref), 'k'); hold on;
% plot(x, real(u), 'r');
plot(x, real(u * 1i * (k)), 'r:');
% plot(real(u))

% plot(k_vec, abs(uu_vec));
% plot(k_vec, angle(uu_vec));

% plot(k_vec2, abs(luu_vec_ref));
% plot(k_vec2, angle(luu_vec_ref));

