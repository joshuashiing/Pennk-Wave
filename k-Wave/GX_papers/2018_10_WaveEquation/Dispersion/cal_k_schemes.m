function cal_k_schemes(f0_m, c0_m, rho, Q, fl, fu, nf, casename)
% Calculate wavenumber given media parameters and frequency ranges
% Author: Guangchi Xing
% Date: 03/03/2018

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

w0 = 2 * pi * f0;
c0 = c0_m * (w0 / w0_m)^gamma;
M0 = rho * c0^2 * (cos(pi*gamma/2))^2;
c = sqrt(M0 / rho);
M = M0 * (1i * w / w0) .^ (2 * gamma);

k = sqrt(rho / M0) * (1i * w / w0) .^ (-gamma) .* w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the wavenumber of approximate wave equation

C_111110 = [-gamma*c*w0; c^2; gamma*c^3/w0; pi*gamma*c; pi*gamma^2*c^2/w0];
C_111100 = [-gamma*c*w0; c^2; gamma*c^3/w0; pi*gamma*c];
C_111010 = [-gamma*c*w0; c^2; gamma*c^3/w0; 10/19*pi*gamma*c^2/w0];
% C_011100 = [(1-17/14*gamma)*c^2; (10/7*gamma+101/37*gamma^2)*c^3/w0; (pi*gamma+39/7*gamma^2)*c];
% C_110100 = [-(3*gamma+23/6*gamma^2)*c*w0; (1+8/3*gamma)*c^2; (pi*gamma+4/3*pi*gamma^2)*c];
C_011100 = [(1)*c^2; (10/7*gamma)*c^3/w0; (pi*gamma)*c];
C_110100 = [-(3*gamma)*c*w0; (1)*c^2; (pi*gamma)*c];

k_111110 = solve_disp(w, k, C_111110, '111110');
k_111100 = solve_disp(w, k, C_111100, '111100');
k_111010 = solve_disp(w, k, C_111010, '111010');
k_011100 = solve_disp(w, k, C_111110, '011100');
k_110100 = solve_disp(w, k, C_111110, '110100');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the resultant wavenumber
save(casename, 'f', 'k', 'k_111110', 'k_111100', 'k_111010', 'k_011100', 'k_110100');