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
f0 = 30;       % Reference frequency for calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the frequency of interest
nf  = 100;
f   = linspace(10, 60, nf);
% f = 40;
w   = 2 * pi * f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the relevant parameters (Kjartansson's model)
gamma = 1 / pi * atan(1 / Q);
M0_m = rho * c0_m^2 * (cos(pi*gamma/2))^2;
w0_m = 2 * pi * f0_m;
c_m = sqrt(M0_m / rho);

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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NCQ Dispersion Relation in Zhu & Harris, 2014
% mod_mech = 'TZ14';
% c = sqrt(M0 / rho);
% % lhs = (w .^ 2) / (c .^ 2);
% % rhs = c0^(2*gamma) * w0^(-2*gamma) * cos(pi*gamma) * k.^(2*gamma+2) + ...
% %     (1i*w) .* c0^(2*gamma-1) * w0^(-2*gamma) * sin(pi*gamma) .* k.^(2*gamma+1);
% 
% k_sol_list = [];
% % nf = 1;
% k_sol_list = zeros(nf);
% for i = 1 : nf
%     i
%     ww = w(i);
%     lhs = ww^2 / c^2; 
%     C1 = c0^(2*gamma) * w0^(-2*gamma) * cos(pi*gamma);
%     C2 = (1i*ww) * c0^(2*gamma-1) * w0^(-2*gamma) * sin(pi*gamma);
%     exp1 = 2 * gamma + 2;
%     exp2 = 2 * gamma + 1;
%     syms kk;
%     k_root = vpasolve(lhs - C1 * kk^exp1 - C2 * kk^exp2 == 0, kk, w0/c0);
%     k_sol(i) = double(k_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Taylor-Fitting (TF17)
% % w^2 = c1 * k + c2 * k^2 + c3 * k^3 
% %       + c4 * (i*w) * k + c5 * (i*w) * k^2 + c6 * (i*w) * k^3
% 
% mod_mech = 'TF17';
% 
% % Generate the B matrix
% B = zeros(6, 6);
% B_tmp = [1; (1 - gamma); (-1/2*gamma + 1/2*gamma^2)];
% B(1:3, 1) = B_tmp * cos(-pi*gamma/2);
% B(4:6, 1) = B_tmp * sin(-pi*gamma/2);
% B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
% B(1:3, 2) = B_tmp * cos(-pi*gamma);
% B(4:6, 2) = B_tmp * sin(-pi*gamma);
% B_tmp = [1; (3 - 3*gamma); (3 - 15/2*gamma + 9/2*gamma^2)];
% B(1:3, 3) = B_tmp * cos(-3/2*pi*gamma);
% B(4:6, 3) = B_tmp * sin(-3/2*pi*gamma);
% B_tmp = [1; (2 - gamma); (1 - 3/2*gamma + 1/2*gamma^2)];
% B(1:3, 4) = B_tmp * cos(pi/2 - pi*gamma/2);
% B(4:6, 4) = B_tmp * sin(pi/2 - pi*gamma/2);
% B_tmp = [1; (3 - 2*gamma); (3 - 5*gamma + 2*gamma^2)];
% B(1:3, 5) = B_tmp * cos(pi/2 - pi * gamma);
% B(4:6, 5) = B_tmp * sin(pi/2 - pi * gamma);
% B_tmp = [1; (4 - 3*gamma); (6 - 21/2*gamma + 9/2*gamma^2)];
% B(1:3, 6) = B_tmp * cos(pi/2 - 3/2*pi*gamma);
% B(4:6, 6) = B_tmp * sin(pi/2 - 3/2*pi*gamma);
% 
% % Solve Linear System
% A = B \ [1; 2; 1; 0; 0; 0];
% 
% mod_mech = 'DT-TF17';
% clear A
% A(1) = -gamma;
% A(2) = 1.0;
% A(3) = gamma;
% A(4) = pi * gamma;
% A(5) = pi * gamma^2;
% A(6) = -3/2 * pi * gamma^4;
% 
% C_k1 = A(1) * c * w0;
% C_k2 = A(2) * c^2;
% C_k3 = A(3) * c^3 / w0;
% C_k4 = A(4) * c;
% C_k5 = A(5) * c^2 / w0;
% C_k6 = A(6) * c^3 / (w0^2);
% 
% % Solve the dispersion relation
% k_sol = zeros(size(w));
% for i = 1 : nf
%     i
%     w_i = w(i);
%     lhs = w_i ^ 2;
%     
%     syms kk;
%     eqn = (lhs == C_k1 * kk + C_k2 * kk^2 + C_k3 * kk^3 + ...
%         C_k4 * (1i*w_i) * kk + C_k5 * (1i*w_i) * kk^2 + C_k6 * (1i*w_i) * kk^3);
%     k_root = vpasolve(eqn, kk, k(i));
%     k_root = double(k_root);
%     [~, i_root] = min(abs(k_root - k(i)));
%     k_sol(i) = k_root(i_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Five-Term (FT17)
% % w^2 = c1 * k + c2 * k^2 + c3 * k^3 
% %       + c4 * (i*w) * k + c5 * (i*w) * k^2
% 
% mod_mech = 'FT17';
% 
% % Generate the B matrix
% B = zeros(6, 5);
% B_tmp = [1; (1 - gamma); (-1/2*gamma + 1/2*gamma^2)];
% B(1:3, 1) = B_tmp * cos(-pi*gamma/2);
% B(4:6, 1) = B_tmp * sin(-pi*gamma/2);
% B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
% B(1:3, 2) = B_tmp * cos(-pi*gamma);
% B(4:6, 2) = B_tmp * sin(-pi*gamma);
% B_tmp = [1; (3 - 3*gamma); (3 - 15/2*gamma + 9/2*gamma^2)];
% B(1:3, 3) = B_tmp * cos(-3/2*pi*gamma);
% B(4:6, 3) = B_tmp * sin(-3/2*pi*gamma);
% B_tmp = [1; (2 - gamma); (1 - 3/2*gamma + 1/2*gamma^2)];
% B(1:3, 4) = B_tmp * cos(pi/2 - pi*gamma/2);
% B(4:6, 4) = B_tmp * sin(pi/2 - pi*gamma/2);
% B_tmp = [1; (3 - 2*gamma); (3 - 5*gamma + 2*gamma^2)];
% B(1:3, 5) = B_tmp * cos(pi/2 - pi * gamma);
% B(4:6, 5) = B_tmp * sin(pi/2 - pi * gamma);
% 
% % Solve Linear System
% A = (B' * B) \ (B' * [1; 2; 1; 0; 0; 0]);
% 
% mod_mech = 'DT-FT17';
% clear A
% A(1) = -gamma;
% A(2) = 1.0;
% A(3) = gamma;
% A(4) = pi * gamma;
% A(5) = pi * gamma^2;
% 
% C_k1 = A(1) * c * w0;
% C_k2 = A(2) * c^2;
% C_k3 = A(3) * c^3 / w0;
% C_k4 = A(4) * c;
% C_k5 = A(5) * c^2 / w0;
% 
% % Solve the dispersion relation
% k_sol = zeros(size(w));
% for i = 1 : nf
%     i
%     w_i = w(i);
%     lhs = w_i ^ 2;
%     
%     syms kk;
%     eqn = (lhs == C_k1 * kk + C_k2 * kk^2 + C_k3 * kk^3 + ...
%         C_k4 * (1i*w_i) * kk + C_k5 * (1i*w_i) * kk^2);
%     k_root = vpasolve(eqn, kk, k(i));
%     k_root = double(k_root);
%     [~, i_root] = min(abs(k_root - k(i)));
%     k_sol(i) = k_root(i_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Three-plus-One (TO17)
% % w^2 = c1 * k + c2 * k^2 + c3 * k^3 + c4 * (i*w) * k^2
% 
% mod_mech = 'TO17';
% 
% % Generate the B matrix
% B = zeros(6, 4);
% B_tmp = [1; (1 - gamma); (-1/2*gamma + 1/2*gamma^2)];
% B(1:3, 1) = B_tmp * cos(-pi*gamma/2);
% B(4:6, 1) = B_tmp * sin(-pi*gamma/2);
% B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
% B(1:3, 2) = B_tmp * cos(-pi*gamma);
% B(4:6, 2) = B_tmp * sin(-pi*gamma);
% B_tmp = [1; (3 - 3*gamma); (3 - 15/2*gamma + 9/2*gamma^2)];
% B(1:3, 3) = B_tmp * cos(-3/2*pi*gamma);
% B(4:6, 3) = B_tmp * sin(-3/2*pi*gamma);
% B_tmp = [1; (3 - 2*gamma); (3 - 5*gamma + 2*gamma^2)];
% B(1:3, 4) = B_tmp * cos(pi/2 - pi * gamma);
% B(4:6, 4) = B_tmp * sin(pi/2 - pi * gamma);
% 
% % Solve Linear System
% A = (B' * B) \ (B' * [1; 2; 1; 0; 0; 0]);
% 
% mod_mech = 'DT-TO17';
% clear A
% A(1) = -gamma;
% A(2) = 1.0;
% A(3) = gamma;
% A(4) = 1.65346981767884 * gamma;
% 
% C_k1 = A(1) * c * w0;
% C_k2 = A(2) * c^2;
% C_k3 = A(3) * c^3 / w0;
% C_k4 = A(4) * c^2 / w0;
% 
% % Solve the dispersion relation
% k_sol = zeros(size(w));
% for i = 1 : nf
%     i
%     w_i = w(i);
%     lhs = w_i ^ 2;
%     
%     syms kk;
%     eqn = (lhs == C_k1 * kk + C_k2 * kk^2 + C_k3 * kk^3 + ...
%         C_k4 * (1i*w_i) * kk^2);
%     k_root = vpasolve(eqn, kk, k(i));
%     k_root = double(k_root);
%     [~, i_root] = min(abs(k_root - k(i)));
%     k_sol(i) = k_root(i_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Three-plus-One (TO18)
% % w^2 = c1 * k + c2 * k^2 + c3 * k^3 + c4 * (i*w) * k
% 
% mod_mech = 'TO18';
% 
% % Generate the B matrix
% B = zeros(6, 4);
% B_tmp = [1; (1 - gamma); (-1/2*gamma + 1/2*gamma^2)];
% B(1:3, 1) = B_tmp * cos(-pi*gamma/2);
% B(4:6, 1) = B_tmp * sin(-pi*gamma/2);
% B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
% B(1:3, 2) = B_tmp * cos(-pi*gamma);
% B(4:6, 2) = B_tmp * sin(-pi*gamma);
% B_tmp = [1; (3 - 3*gamma); (3 - 15/2*gamma + 9/2*gamma^2)];
% B(1:3, 3) = B_tmp * cos(-3/2*pi*gamma);
% B(4:6, 3) = B_tmp * sin(-3/2*pi*gamma);
% B_tmp = [1; (2 - gamma); (1 - 3/2*gamma + 1/2*gamma^2)];
% B(1:3, 4) = B_tmp * cos(pi/2 - pi*gamma/2);
% B(4:6, 4) = B_tmp * sin(pi/2 - pi*gamma/2);
% 
% % Solve Linear System
% A = (B' * B) \ (B' * [1; 2; 1; 0; 0; 0]);
% 
% mod_mech = 'DT-TO18';
% clear A
% A(1) = -gamma;
% A(2) = 1.0;
% A(3) = gamma;
% A(4) = pi * gamma;
% 
% C_k1 = A(1) * c * w0;
% C_k2 = A(2) * c^2;
% C_k3 = A(3) * c^3 / w0;
% C_k4 = A(4) * c;
% 
% % Solve the dispersion relation
% k_sol = zeros(size(w));
% for i = 1 : nf
%     i
%     w_i = w(i);
%     lhs = w_i ^ 2;
%     
%     syms kk;
%     eqn = (lhs == C_k1 * kk + C_k2 * kk^2 + C_k3 * kk^3 + ...
%         C_k4 * (1i*w_i) * kk);
%     k_root = vpasolve(eqn, kk, k(i));
%     k_root = double(k_root);
%     [~, i_root] = min(abs(k_root - k(i)));
%     k_sol(i) = k_root(i_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Taylor-Fitting (MO18)
% % w^2 = c1 * k^2 + c2 * k^3 + c3 * (i*w) * k
% 
% mod_mech = 'MO18';
% 
% % Generate the B matrix
% B = zeros(6, 3);
% B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
% B(1:3, 1) = B_tmp * cos(-pi*gamma);
% B(4:6, 1) = B_tmp * sin(-pi*gamma);
% B_tmp = [1; (3 - 3*gamma); (3 - 15/2*gamma + 9/2*gamma^2)];
% B(1:3, 2) = B_tmp * cos(-3/2*pi*gamma);
% B(4:6, 2) = B_tmp * sin(-3/2*pi*gamma);
% B_tmp = [1; (2 - gamma); (1 - 3/2*gamma + 1/2*gamma^2)];
% B(1:3, 3) = B_tmp * cos(pi/2 - pi*gamma/2);
% B(4:6, 3) = B_tmp * sin(pi/2 - pi*gamma/2);
% 
% % Solve Linear System
% A = (B' * B) \ (B' * [1; 2; 1; 0; 0; 0]);
% 
% % mod_mech = 'DT-MO18';
% % clear A
% % A(1) = 1.0;
% % A(2) = 10/7 * gamma;
% % A(3) = pi * gamma;
% 
% mod_mech = 'DT2-MO18';
% clear A
% A(1) = 1.0 - 17/14 * gamma;
% A(2) = 10/7 * gamma + 101/37 * gamma^2;
% A(3) = pi * gamma + 39/7 * gamma^2;
% 
% 
% C_k1 = A(1) * c^2;
% C_k2 = A(2) * c^3 / w0;
% C_k3 = A(3) * c;
% 
% % Solve the dispersion relation
% k_sol = zeros(size(w));
% for i = 1 : nf
%     i
%     w_i = w(i);
%     lhs = w_i ^ 2;
%     
%     syms kk;
%     eqn = (lhs == C_k1 * kk^2 + C_k2 * kk^3 + C_k3 * (1i*w_i) * kk);
%     k_root = vpasolve(eqn, kk, k(i));
%     k_root = double(k_root);
%     [~, i_root] = min(abs(k_root - k(i)));
%     k_sol(i) = k_root(i_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Taylor-Fitting (MT18)
% % w^2 = c1 * k + c2 * k^2 + c3 * (i*w) * k
% 
% mod_mech = 'MT18';
% 
% % Generate the B matrix
% B = zeros(6, 3);
% B_tmp = [1; (1 - gamma); (-1/2*gamma + 1/2*gamma^2)];
% B(1:3, 1) = B_tmp * cos(-pi*gamma/2);
% B(4:6, 1) = B_tmp * sin(-pi*gamma/2);
% B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
% B(1:3, 2) = B_tmp * cos(-pi*gamma);
% B(4:6, 2) = B_tmp * sin(-pi*gamma);
% B_tmp = [1; (2 - gamma); (1 - 3/2*gamma + 1/2*gamma^2)];
% B(1:3, 3) = B_tmp * cos(pi/2 - pi*gamma/2);
% B(4:6, 3) = B_tmp * sin(pi/2 - pi*gamma/2);
% 
% % Solve Linear System
% A = (B' * B) \ (B' * [1; 2; 1; 0; 0; 0]);
% 
% % mod_mech = 'DT-MT18';
% % clear A
% % A(1) = -3 * gamma;
% % A(2) = 1.0;
% % A(3) = pi * gamma;
% 
% mod_mech = 'DT2-MT18';
% clear A
% A(1) = -3 * gamma - 23/6 * gamma^2;
% A(2) = 1.0 + 8/3 * gamma;
% A(3) = pi * gamma + 4/3*pi * gamma^2;
% 
% C_k1 = A(1) * c * w0;
% C_k2 = A(2) * c^2;
% C_k3 = A(3) * c;
% 
% % Solve the dispersion relation
% k_sol = zeros(size(w));
% for i = 1 : nf
%     i
%     w_i = w(i);
%     lhs = w_i ^ 2;
%     
%     syms kk;
%     eqn = (lhs == C_k1 * kk + C_k2 * kk^2 + C_k3 * (1i*w_i) * kk);
%     k_root = vpasolve(eqn, kk, k(i));
%     k_root = double(k_root);
%     [~, i_root] = min(abs(k_root - k(i)));
%     k_sol(i) = k_root(i_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taylor-Fitting (FN18)
% w^2 = c1 * k^2 + c2 * k^3 + c3 * (i*w) * k

mod_mech = 'FN18';

% Generate the B matrix
B = zeros(6, 2);

B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
B(1:3, 1) = B_tmp * cos(-pi*gamma);
B(4:6, 1) = B_tmp * sin(-pi*gamma);
B_tmp = [1; (2 - gamma); (1 - 3/2*gamma + 1/2*gamma^2)];
B(1:3, 2) = B_tmp * cos(pi/2 - pi*gamma/2);
B(4:6, 2) = B_tmp * sin(pi/2 - pi*gamma/2);


% Solve Linear System
A = (B' * B) \ (B' * [1; 2; 1; 0; 0; 0]);

% mod_mech = 'DT-FN18';
% clear A
% A(1) = 1.0;
% A(2) = 10/7 * gamma;
% A(3) = pi * gamma;

% mod_mech = 'DT2-FN18';
% clear A
% A(1) = 1.0 - 17/14 * gamma;
% A(2) = 10/7 * gamma + 101/37 * gamma^2;
% A(3) = pi * gamma + 39/7 * gamma^2;


C_k1 = A(1) * c^2;
C_k2 = A(2) * c;

% Solve the dispersion relation
k_sol = zeros(size(w));
for i = 1 : nf
    i
    w_i = w(i);
    lhs = w_i ^ 2;
    
    syms kk;
    eqn = (lhs == C_k1 * kk^2 + C_k2 * (1i*w_i) * kk);
    k_root = vpasolve(eqn, kk, k(i));
    k_root = double(k_root);
    [~, i_root] = min(abs(k_root - k(i)));
    k_sol(i) = k_root(i_root);
end
cp_sol = w ./ real(k_sol);
alpha_sol = -imag(k_sol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% 
% % Tes t plot for complex value
% val = k_m;
% plot(f, real(val), 'k', 'linewidth', 2); hold on;
% plot(f, imag(val), 'r', 'linewidth', 2);
% legend('Real', 'Imaginary');

% % Test plot for real value
% val = alpha_m;
% plot(f, val, 'k', 'linewidth', 2);

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

