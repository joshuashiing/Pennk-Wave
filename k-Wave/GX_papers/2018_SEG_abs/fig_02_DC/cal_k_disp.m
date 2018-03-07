function cal_k_disp(f0_m, c0_m, rho, Q, fl, fu, nf, casename)
% Calculate wavenumber given media parameters and frequency ranges
% Author: Guangchi Xing
% Date: 03/03/2018

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set up the parameters for Pierre shale (Wuenschel, 1965)
% f0_m = 200;     % Reference frequency of the medium (Hz)
% c0_m = 2089;    % Reference phase velocity (m/s)
% rho  = 2200;    % Density (kg/m^3)
% Q  = 32.5;      % Quality factor
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set up the frequency of interest
% fl = 10;        % Lower bound of frequency range (Hz)
% fu = 100;       % Upper bound of frequency range (Hz)
% nf = 100;       % Number of frequency sample points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the frequency axis
f0 = (fl + fu) / 2;         % Reference frequency
f = linspace(fl, fu, nf);   % frequency axis
w = 2 * pi * f;             % angular frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the relevant parameters (Kjartansson's model)
gamma = 1 / pi * atan(1 / Q);
M0_m = rho * c0_m^2 * (cos(pi*gamma/2))^2;
w0_m = 2 * pi * f0_m;
% c_m = sqrt(M0_m / rho);

w0 = 2 * pi * f0;
c0 = c0_m * (w0 / w0_m)^gamma;
M0 = rho * c0^2 * (cos(pi*gamma/2))^2;
c = sqrt(M0 / rho);
% alpha0 = 1 / c0 * w0^gamma * tan(pi*gamma/2);
M = M0 * (1i * w / w0) .^ (2 * gamma);

k = sqrt(rho / M0) * (1i * w / w0) .^ (-gamma) .* w;
k_r = real(k);
k_i = imag(k);
vc = sqrt(M0 / rho) * (1i * w / w0) .^ gamma;
cp = c0 * (w / w0) .^ gamma;
alpha = w ./ cp * tan(pi * gamma / 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate wave equation coefficients

% Normalization coeffecient
C_norm = [c*w0; c^2; c^3/w0; c; c^2/w0; c^3/(w0^2)];

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

% Solve for A vector and obtain C vector
C_111111 = solve_coef(B, C_norm, [1, 2, 3, 4, 5, 6]);
C_111110 = solve_coef(B, C_norm, [1, 2, 3, 4, 5]);
C_111100 = solve_coef(B, C_norm, [1, 2, 3, 4]);
C_111010 = solve_coef(B, C_norm, [1, 2, 3, 5]);
C_011100 = solve_coef(B, C_norm, [2, 3, 4]);
C_110100 = solve_coef(B, C_norm, [1, 2, 4]);

% Closed form C vector obtained by Taylor expansion
C_c_111111 = [-gamma*c*w0; c^2; gamma*c^3/w0; pi*gamma*c; pi*gamma^2*c^2/w0; -3/2*pi*gamma^4*c^3/w0^2];
C_c_111110 = [-gamma*c*w0; c^2; gamma*c^3/w0; pi*gamma*c; pi*gamma^2*c^2/w0];
C_c_111100 = [-gamma*c*w0; c^2; gamma*c^3/w0; pi*gamma*c];
C_c_011100 = [(1-17/14*gamma)*c^2; (10/7*gamma+101/37*gamma^2)*c^3/w0; (pi*gamma+39/7*gamma^2)*c];
C_c_011100_b = [c^2; 10/7*gamma*c^3/w0; pi*gamma*c];
C_c_110100 = [-(3*gamma+23/6*gamma^2)*c*w0; (1+8/3*gamma)*c^2; (pi*gamma+4/3*pi*gamma^2)*c];
C_c_110100_b = [-3*gamma*c*w0; c^2; pi*gamma*c];

% Numerically solve dispersion relation for wavenumber
k_111111 = solve_disp(w, k, C_111111, '111111');
k_111110 = solve_disp(w, k, C_111110, '111110');
k_111100 = solve_disp(w, k, C_111100, '111100');
k_111010 = solve_disp(w, k, C_111010, '111010');
k_011100 = solve_disp(w, k, C_011100, '011100');
k_110100 = solve_disp(w, k, C_110100, '110100');

k_c_111111 = solve_disp(w, k, C_c_111111, '111111');
k_c_111110 = solve_disp(w, k, C_c_111110, '111110');
k_c_111100 = solve_disp(w, k, C_c_111100, '111100');
k_c_011100 = solve_disp(w, k, C_c_011100, '011100');
k_c_011100_b = solve_disp(w, k, C_c_011100_b, '011100');
k_c_110100 = solve_disp(w, k, C_c_110100, '110100');
k_c_110100_b = solve_disp(w, k, C_c_110100_b, '110100');

% Save the resultant wavenumber
save([casename, '_k.mat'], 'f', 'k', ...
    'k_111111', 'k_111110', 'k_111100', 'k_111010', 'k_011100', 'k_110100', ...
    'k_c_111111', 'k_c_111110', 'k_c_111100', 'k_c_011100', 'k_c_011100_b', 'k_c_110100', 'k_c_110100_b');
