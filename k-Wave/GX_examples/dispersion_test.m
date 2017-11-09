clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the characteristic parameters of exploration seismology
% c0  |  2  - 4 km/s
% f   |  15 - 40 Hz
% rho |  2  - 2.5 g/cm^3
% Q   |  20 - 100

f0_m = 200;     % Reference frequency of the medium
c0_m = 2089;    % Reference phase velocity
rho = 2200;
Q   = 32;
f0 = 25;       % Reference frequency for calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the frequency of interest
nf  = 100;
f   = linspace(10, 60, nf);
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



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Test for the Kramers-Kronig Relation for 0 < y < 1
% lhs = 1 ./ cp;
% rhs = alpha0 * tan(y * pi / 2) * w .^ (y - 1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NCQ Dispersion Relation in Zhu & Harris, 2014
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NCQ Dispersion for Zhu 2017 Draft
% a0 = pi * gamma * w0^gamma / (2 * c0);
% % a0 = pi * gamma / (2 * c0 * w0^gamma);
% % 
% % % Eqn (6)
% % lhs = k;
% % rhs = w ./ cp - 1i * a0 * w;
% 
% % % Eqn (9)
% % rhs = (1+gamma)^2 / (c0^2) * w.^2 - ...
% %     2 * gamma * (1+gamma) * c0 / w0 * w.^3 / (c0^3) - ...
% %     2 * a0 * w / cp .* (1i*w);
% % k_sol = sqrt(rhs);
% % cp_sol = w ./ real(k_sol);
% % alpha_sol = -imag(k_sol);
% 
% % Eqn (10)
% k_sol = zeros(size(w));
% for i = 1 : nf
%     i
%     w_i = w(i);
%     lhs = - (1 + gamma)^2 / (c0^2) * w_i^2;
%     c_k1 = -2 * a0 * (1i * w_i);
%     c_k2 = -1;
%     c_k3 = -2 * gamma * (1+gamma) * c0 / w0;
%     
%     syms kk;
%     eqn = (lhs == c_k1*kk + c_k2*kk^2 + c_k3*kk^3);
%     k_root = vpasolve(eqn, kk, k(i));
%     k_root = double(k_root);
%     [~, i_root] = min(abs(k_root - k(i)));
%     k_sol(i) = k_root(i_root);
% end
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);
% 
% % % Eqn (10) Coefficient
% % C_k2 = c0^2 / (1+gamma)^2;
% % C_k3 = 2*gamma*(1+gamma)*c0 / w0 * c0^2 / (1+gamma)^2;
% % C_k4 = 2 * a0;
% % C_k_Zhu = [0; C_k2; C_k3; C_k4; 0; 0];
% % clear C_k2 C_k3 C_k4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My Dispersion
% c = sqrt(M0 / rho);
% % lhs = w.^2 / (c^2);
% % rhs = (1i)^(2*gamma) * k.^2 .* (gamma*(2*gamma-1)*(w/w0).^2 + 4*gamma*(1-gamma)*(w/w0) + (gamma-1)*(2*gamma-1));
% lhs = w / c;
% rhs = (1i)^gamma * k .* (sqrt(gamma*(2*gamma-1)) * (w/w0) + 2*gamma*(1-gamma)/sqrt(gamma*(2*gamma-1)));
% lhs = w.^2;
% rhs = (M0 / rho)^(1/(1-gamma)) * exp(1i*pi*gamma/(1-gamma)) * w0^(-2*gamma/(1-gamma)) * k.^(2/(1-gamma));
% lhs = k .^ (2 / (1-gamma));
% rhs = (w0/c0)^(-2/(gamma-1)) + 2/(1-gamma) * (w0/c0)^((1+gamma)/1-gamma) * (k - w0/c0) + (1+gamma)/((1-gamma)^2) * (w0/c0)^(-2*gamma/(gamma-1)) * (k-w0/c0).^2;
% rhs = (w0/c0)^(-2/(gamma-1)) + 2/(1-gamma) * (w0/c0)^((1+gamma)/1-gamma) * (k - w0/c0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taylor Fitting Dispersion

% Generate the B matrix
B = zeros(6, 6);
B_tmp = [1; (1 - gamma); (-1/2*gamma + 1/2*gamma^2)];
B(1:3, 1) = B_tmp * cos(-pi*gamma/2);
B(4:6, 1) = B_tmp * sin(-pi*gamma/2);
B_tmp = [1; (2 - 2*gamma); (1 - 3*gamma + 2*gamma^2)];
B(1:3, 2) = B_tmp * cos(-pi*gamma);
B(4:6, 2) = B_tmp * sin(-pi*gamma);
B_tmp = [1; (3 - 3*gamma); (3 - 15/2*gamma + 9/2*gamma^2)];
B(1:3, 3) = B_tmp * cos(-3/2*pi*gamma);
B(4:6, 3) = B_tmp * sin(-3/2*pi*gamma);
B_tmp = [1; (2 - gamma); (1 - 3/2*gamma + 1/2*gamma^2)];
B(1:3, 4) = B_tmp * cos(pi/2 - pi*gamma/2);
B(4:6, 4) = B_tmp * sin(pi/2 - pi*gamma/2);
B_tmp = [1; (3 - 2*gamma); (3 - 5*gamma + 2*gamma^2)];
B(1:3, 5) = B_tmp * cos(pi/2 - pi * gamma);
B(4:6, 5) = B_tmp * sin(pi/2 - pi * gamma);
B_tmp = [1; (4 - 3*gamma); (6 - 21/2*gamma + 9/2*gamma^2)];
B(1:3, 6) = B_tmp * cos(pi/2 - 3/2*pi*gamma);
B(4:6, 6) = B_tmp * sin(pi/2 - 3/2*pi*gamma);

% Check Taylor expansion
x = (w - w0) ./ w0;
x_vec = [ones(size(x)); x; x.^2];
k1_te = (B(1:3, 1)' * x_vec + 1i * B(4:6, 1)' * x_vec) * w0 / c;
k2_te = (B(1:3, 2)' * x_vec + 1i * B(4:6, 2)' * x_vec) * w0^2 / (c^2);
k3_te = (B(1:3, 3)' * x_vec + 1i * B(4:6, 3)' * x_vec) * w0^3 / (c^3);
k4_te = (B(1:3, 4)' * x_vec + 1i * B(4:6, 4)' * x_vec) * w0^2 / c;
k5_te = (B(1:3, 5)' * x_vec + 1i * B(4:6, 5)' * x_vec) * w0^3 / (c^2);
k6_te = (B(1:3, 6)' * x_vec + 1i * B(4:6, 6)' * x_vec) * w0^4 / (c^3);


% Solve B * A = [1; 2; 1; 0; 0; 0] and assign k coefficients
A = B \ [1; 2; 1; 0; 0; 0];
% A(6) = 0;
% A(5) = 0;
% A(1) = 0;
C_k1 = A(1) * c * w0;
C_k2 = A(2) * c^2;
C_k3 = A(3) * c^3 / w0;
C_k4 = A(4) * c;
C_k5 = A(5) * c^2 / w0;
C_k6 = A(6) * c^3 / (w0^2);

C_k = [C_k1; C_k2; C_k3; C_k4; C_k5; C_k6];

% Solve the dispersion relation
k_sol = zeros(size(w));
for i = 1 : nf
    i
    w_i = w(i);
    lhs = w_i ^ 2;
    
    syms kk;
    eqn = (lhs == C_k1 * kk + C_k2 * kk^2 + C_k3 * kk^3 + ...
        C_k4 * (1i*w_i) * kk + C_k5 * (1i*w_i) * kk^2 + C_k6 * (1i*w_i) * kk^3);
    k_root = vpasolve(eqn, kk, k(i));
    k_root = double(k_root);
    [~, i_root] = min(abs(k_root - k(i)));
    k_sol(i) = k_root(i_root);
end
cp_sol = w ./ real(k_sol);
alpha_sol = -imag(k_sol);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Play with (w/w0)^(2*gamma)
% lhs = (w ./ w0) .^ (2*gamma);
% dww0 = (w - w0) ./ w0;
% % rhs = 1 + 2*gamma*dww0 +...
% %     1/2 * (2*gamma) * (2*gamma-1) * dww0.^2 +...
% %     1/6 * (2*gamma) * (2*gamma-1) * (2*gamma-2) * dww0.^3;
% rhs = 1 + 2*gamma*dww0 +...
%     1/2 * (2*gamma) * (2*gamma-1) * dww0.^2;
% 
% % Modified wave equation
% c = sqrt(M0 / rho);
% chunk = gamma * (2*gamma-1) / (w0^2) * w.^2 + ...
%     4 * gamma * (1-gamma) / w0 * w + (gamma-1) * (2*gamma-1);
% k_sol = (1i)^(-gamma) * w / c ./ sqrt(chunk);
% cp_sol = w ./ real(k_sol);
% alpha_sol = -imag(k_sol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

% % Tes t plot for complex value
% val = k_m;
% plot(f, real(val), 'k', 'linewidth', 2); hold on;
% plot(f, imag(val), 'r', 'linewidth', 2);
% legend('Real', 'Imaginary');

% % Test plot for real value
% val = alpha_m;
% plot(f, val, 'k', 'linewidth', 2);

% Test for two values
val1 = alpha;
val2 = alpha_sol;
subplot(211);
plot(f, val1, 'k', 'linewidth', 2); hold on;
plot(f, val2, 'r--', 'linewidth', 2);
subplot(212);
plot(f, (val1 - val2) ./ val1 * 100, 'r', 'linewidth', 2);


