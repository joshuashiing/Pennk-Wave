clear;
% clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the characteristic parameters of exploration seismology
% c0  |  2  - 4 km/s
% f   |  15 - 40 Hz
% rho |  2  - 2.5 g/cm^3
% Q   |  20 - 100

f0_m = 200;     % Reference frequency of the medium
c0_m = 2089;    % Reference phase velocity
rho = 2200;
Q   = 32.5;
f0 = 300;       % Reference frequency for calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the frequency of interest
nf  = 100;
f   = linspace(10, 600, nf);
% f = 40;
w   = 2 * pi * f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the relevant parameters (Kjartansson's model)
gamma = 1 / pi * atan(1 / Q);
M0_m = rho * c0_m^2 * (cos(pi*gamma/2))^2;
w0_m = 2 * pi * f0_m;
c_m = sqrt(M0_m / rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapt to the acoustic community:
% alpha = alpha0 * w^y
alpha0_m = 1 / c0_m * w0_m ^ gamma * tan(pi * gamma / 2);
y = 1 - gamma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the dispersion relations
M_m = M0_m * (1i * w / w0_m) .^ (2 * gamma);
k_m = sqrt(rho / M0_m) * (1i * w / w0_m) .^ (-gamma) .* w;    % Complex wavenum.
k_m_r = real(k_m);
k_m_i = imag(k_m);
vc_m = sqrt(M0_m / rho) * (1i * w / w0_m) .^ gamma;           % Complex velocity
cp_m = c0_m * (w / w0_m) .^ gamma;                            % Phase velocity
alpha_m = w ./ cp_m * tan(pi * gamma / 2);                  % Attenuation Fac.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a specific reference frequency f0, calculate the relevant
% parameters for Kjartansson's model
w0 = 2 * pi * f0;
c0 = c0_m * (w0 / w0_m)^gamma;
M0 = rho * c0^2 * (cos(pi*gamma/2))^2;
c = sqrt(M0 / rho);
alpha0 = 1 / c0 * w0^gamma * tan(pi*gamma/2);
M = M0 * (1i * w / w0) .^ (2 * gamma);
k = sqrt(rho / M0) * (1i * w / w0) .^ (-gamma) .* w;
k_r = real(k);
k_i = imag(k);
vc = sqrt(M0 / rho) * (1i * w / w0) .^ gamma;
cp = c0 * (w / w0) .^ gamma;
alpha = w ./ cp * tan(pi * gamma / 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDTF(Finite-Difference Taylor-Fitting) 8-th order

mod_mech = 'PF17';

m_max = 4;
n_max = m_max;
n_x_n = 4;

x = (w - w0) ./ w0;
x_n_vec = transpose(linspace(-1, 1, n_x_n + 1));
x_n_vec = x_n_vec(2:end);
x_n_vec = transpose(linspace(10, 600, n_x_n)) / f0 - 1;
% x_n_vec = [-0.5; 1];

m_vec = 1 : m_max;
n_vec = 1 : n_max;

B = zeros(n_x_n * 2, m_max + n_max);

% for i_x_n = 1 : n_x_n
%     x_n = x_n_vec(i_x_n);
%     B(i_x_n, 1 : m_max) = ...
%         w0.^(2 * m_vec * gamma) .* cos(-m_vec * pi * gamma) .* (1 + x_n).^(2*m_vec*(1-gamma));
%     B(n_x_n + i_x_n, 1 : m_max) = ...
%         w0.^(2 * m_vec * gamma) .* sin(-m_vec * pi * gamma) .* (1 + x_n).^(2*m_vec*(1-gamma));
%     B(i_x_n, (m_max+1) : (m_max+n_max)) = ...
%         w0.^(2 * n_vec * gamma) .* cos(pi/2 - n_vec * pi * gamma) .* (1 + x_n).^(2*n_vec*(1-gamma)+1);
%     B(n_x_n + i_x_n, (m_max+1) : (m_max+n_max)) = ...
%         w0.^(2 * n_vec * gamma) .* sin(pi/2 - n_vec * pi * gamma) .* (1 + x_n).^(2*n_vec*(1-gamma)+1);
% end

for i_x_n = 1 : n_x_n
    x_n = x_n_vec(i_x_n);
    B(i_x_n, 1 : m_max) = ...
        cos(-m_vec * pi * gamma) .* (1 + x_n).^(2*m_vec*(1-gamma));
    B(n_x_n + i_x_n, 1 : m_max) = ...
        sin(-m_vec * pi * gamma) .* (1 + x_n).^(2*m_vec*(1-gamma));
    B(i_x_n, (m_max+1) : (m_max+n_max)) = ...
        cos(pi/2 - n_vec * pi * gamma) .* (1 + x_n).^(2*n_vec*(1-gamma)+1);
    B(n_x_n + i_x_n, (m_max+1) : (m_max+n_max)) = ...
        sin(pi/2 - n_vec * pi * gamma) .* (1 + x_n).^(2*n_vec*(1-gamma)+1);
end

rhs_node = [((1 + x_n_vec).^2); zeros(size(x_n_vec))];
% cd_vec = (B' * B) \ (B' * rhs_node);
cd_vec = B \ rhs_node;

c_vec0 = cd_vec(1 : m_max);
d_vec0 = cd_vec((m_max+1) : end);

lhs_x = (1 + x).^2;
rhs_x = zeros(size(x));

for m = 1 : m_max
    cm = c_vec0(m);
    rhs_x = rhs_x + cm * (1+x).^(2*m*(1-gamma)) * exp(1i*(-m*pi*gamma));
end
for n = 1 : n_max
    dn = d_vec0(n);
    rhs_x = rhs_x + dn * (1+x).^(2*n*(1-gamma)+1) * exp(1i*(pi/2-n*pi*gamma));
end

c_vec = c_vec0 .* w0.^(2-2*m_vec') .* c.^(2*m_vec');
d_vec = d_vec0 .* w0.^(1-2*n_vec') .* c.^(2*n_vec');

c_k1 = c_vec(1);
d_k1 = d_vec(1);
c_k2 = c_vec(2);
d_k2 = d_vec(2);
c_k3 = c_vec(3);
d_k3 = d_vec(3);
c_k4 = c_vec(4);
d_k4 = d_vec(4);

% w2_rhs = zeros(4, length(x));
% w2_rhs(1, :) = c_k1 * k.^2;
% w2_rhs(2, :) = d_k1 * (1i*w) .* k.^2;
% w2_rhs(3, :) = c_k2 * k.^4;
% w2_rhs(4, :) = d_k2 * (1i*w) .* k.^4;
% w2_rhs(5, :) = c_k3 * k.^6;
% w2_rhs(6, :) = d_k3 * (1i*w) .* k.^6;
% w2_rhs = sum(w2_rhs);


% Solve the dispersion relation
k_sol = zeros(size(w));

% c_k1 = c_vec(1);
% c_k2 = c_vec(2);
% d_k1 = d_vec(1);
% d_k2 = d_vec(2);

for i = 1 : nf
    i
    w_i = w(i);
    lhs = w_i ^ 2;
    
%     % test
%     k_i = k(i);
%     lhs = c_k1 * k_i^2 + c_k2 * k_i^4 + ...
%         d_k1 * (1i*w_i) * k_i^2 + d_k2 * (1i*w_i) * k_i^4;
%     % test
    
    syms kk;
    eqn = (lhs == c_k1 * kk^2 + c_k2 * kk^4 + c_k3 * kk^6 + c_k4 * kk^8 + ...
        d_k1 * (1i*w_i) * kk^2 + d_k2 * (1i*w_i) * kk^4 + d_k3 * (1i*w_i) * kk^6) + d_k4 * (1i*w_i) * kk^8;
%     eqn = (lhs == c_k1 * kk^2 + c_k2 * kk^4 + c_k3 * kk^6 + ...
%         d_k1 * (1i*w_i) * kk^2 + d_k2 * (1i*w_i) * kk^4 + d_k3 * (1i*w_i) * kk^6);
%     eqn = (lhs == c_k1 * kk^2 + c_k2 * kk^4 + ...
%         d_k1 * (1i*w_i) * kk^2 + d_k2 * (1i*w_i) * kk^4);
%     eqn = (lhs == C_k1 * kk^2 + C_k2 * kk^4 + C_k3 * kk^6 + + C_k4 * kk^8 + ...
%         C_k5 * (1i*w_i) * kk^2 + C_k6 * (1i*w_i) * kk^4 + C_k7 * (1i*w_i) * kk^6 + C_k8 * (1i*w_i) * kk^8);
%     eqn = (lhs == C_k1 * kk + C_k2 * kk^2 + C_k3 * kk^3 + ...
%         C_k4 * (1i*w_i) * kk + C_k5 * (1i*w_i) * kk^2 + C_k3 * (1i*w_i) * kk^3);
    k_root = vpasolve(eqn, kk, k(i));
    k_root = double(k_root);
    [~, i_root] = min(abs(k_root - k(i)));
%     [~, i_root] = min(abs(k_root - k(i).^2));
    k_sol(i) = k_root(i_root);
%     k_sol(i) = sqrt(k_root(i_root));
end
cp_sol = w ./ real(k_sol);
alpha_sol = -imag(k_sol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;

% % Test for two values
% val1 = w.^2;
% val2 = abs(w2_rhs);
% subplot(211);
% plot(f, val1, 'k', 'linewidth', 2); hold on;
% plot(f, val2, 'r--', 'linewidth', 2);
% subplot(212);
% plot(f, (val2 - val1) ./ val1 * 100, 'r', 'linewidth', 2);



% Test for two groups of real values
val_a1 = cp;
val_a2 = cp_sol;
val_b1 = alpha;
val_b2 = alpha_sol;
subplot(221);
plot(f, val_a1, 'k', 'linewidth', 2); hold on;
plot(f, val_a2, 'r--', 'linewidth', 2);
ylabel('Cp (m/s)')
set(gca, 'fontsize', 14);
title(['Q = ', num2str(Q), ' | f0 = ', num2str(f0), 'Hz'], 'fontsize', 14); 

subplot(223);
plot(f, (val_a2 - val_a1) ./ val_a1 * 100, 'r', 'linewidth', 2);
ylabel('Cp Discrep. (%)')
xlabel('Frequency (Hz)');
set(gca, 'fontsize', 14);

subplot(222);
plot(f, val_b1, 'k', 'linewidth', 2); hold on;
plot(f, val_b2, 'r--', 'linewidth', 2);
ylabel('Alpha (m/s)')
set(gca, 'fontsize', 14);
title(['Red Dashed Line: ', mod_mech], 'fontsize', 14); 

subplot(224);
plot(f, (val_b2 - val_b1) ./ val_b1 * 100, 'r', 'linewidth', 2);
ylabel('Alpha Discrep. (%)')
xlabel('Frequency (Hz)');
set(gca, 'fontsize', 14);

