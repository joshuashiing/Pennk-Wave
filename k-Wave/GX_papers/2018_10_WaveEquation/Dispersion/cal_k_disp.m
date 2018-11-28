function cal_k_disp(f0_m, c0_m, rho, Q, fl, fu, nf, casename)
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
k_111110 = solve_disp(w, k, C_111110, '111110');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the resultant wavenumber
save(casename, 'f', 'k', 'k_111110');