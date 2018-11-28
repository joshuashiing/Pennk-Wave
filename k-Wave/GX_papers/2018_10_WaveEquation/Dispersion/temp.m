cc = C_111100;
kk = k_111100;
rhs = cc(1) * kk + cc(2) * kk.^2 + cc(3) * kk.^3 + cc(4) * kk .* (1i * w);
lhs = w.^2;

cc = C_111010;
kk = k_111010;
rhs = cc(1) * kk + cc(2) * kk.^2 + cc(3) * kk.^3 + cc(4) * kk.^2 .* (1i * w);
lhs = w.^2;

cc = C_011100;
kk = k_011100;
rhs = cc(1) * kk.^2 + cc(2) * kk.^3 + cc(3) * kk .* (1i * w);
lhs = w.^2;







subplot(411);
plot(lhs, 'k', 'linewidth', 2); hold on;
plot(abs(rhs), 'r--', 'linewidth', 2);

subplot(412);
plot(lhs - abs(rhs), 'k', 'linewidth', 2); hold on;

subplot(413);
plot(f, w ./ real(k), 'k', 'linewidth', 2); hold on;
plot(f, w ./ real(kk), 'r--', 'linewidth', 2);
subplot(414);

plot(f, -imag(k), 'k', 'linewidth', 2); hold on;
plot(f, -imag(kk), 'r--', 'linewidth', 2);