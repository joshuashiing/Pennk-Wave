% Define different cases (media parameters and frequency range) and 
% calculate the wavenumber dispersion using code cal_DC.m
% Then Plot by calling plot_DC.m
% Author: Guangchi Xing
% Date: 03/03/2018

clear;
clc

plot_DC('case01', '', 1, 'k', 2);
plot_DC('case01', '111111', 1, 'b:', 2);
plot_DC('case01', 'c_111110', 1, 'r--', 2);
plot_DC('case02', '', 1, 'k', 2);
plot_DC('case02', '111111', 1, 'b:', 2);
plot_DC('case02', 'c_111110', 1, 'r--', 2);
plot_DC('case03', '', 1, 'k', 2);
plot_DC('case03', '111111', 1, 'b:', 2);
plot_DC('case03', 'c_111110', 1, 'r--', 2);
figure(1);
subplot(221); legend('Theoretical', '111111', '111110 close');
subplot(223); legend('111111', '111110 close');
subplot(222); legend('Theoretical', '111111', '111110 close');
subplot(224); legend('111111', '111110 close');


plot_DC('case04', '', 2, 'k', 2);
plot_DC('case04', '111111', 2, 'b:', 2);
plot_DC('case04', 'c_111110', 2, 'r--', 2);
plot_DC('case05', '', 2, 'k', 2);
plot_DC('case05', '111111', 2, 'b:', 2);
plot_DC('case05', 'c_111110', 2, 'r--', 2);
plot_DC('case06', '', 2, 'k', 2);
plot_DC('case06', '111111', 2, 'b:', 2);
plot_DC('case06', 'c_111110', 2, 'r--', 2);
figure(2);
subplot(221); legend('Theoretical', '111111', '111110 close');
subplot(223); legend('111111', '111110 close');
subplot(222); legend('Theoretical', '111111', '111110 close');
subplot(224); legend('111111', '111110 close');


plot_DC('case01', '', 3, 'k', 2);
plot_DC('case01', '111110', 3, 'b:', 2);
plot_DC('case01', '111100', 3, 'r:', 2);
plot_DC('case01', '111010', 3, 'g:', 2);
plot_DC('case01', '011100', 3, 'c:', 2);
plot_DC('case01', '110100', 3, 'm:', 2);
figure(3);
subplot(221); legend('Theoretical', '111110', '111100', '111010', '011100', '110100');
subplot(223); legend('111110', '111100', '111010', '011100', '110100');
subplot(222); legend('Theoretical', '111110', '111100', '111010', '011100', '110100');
subplot(224); legend('111110', '111100', '111010', '011100', '110100');


plot_DC('case01', '', 4, 'k', 2);
plot_DC('case01', 'c_111110', 4, 'b--', 2);
plot_DC('case01', 'c_111100', 4, 'r--', 2);
plot_DC('case01', 'c_011100', 4, 'c--', 2);
plot_DC('case01', 'c_011100_b', 4, 'y--', 2);
plot_DC('case01', 'c_110100', 4, 'm--', 2);
plot_DC('case01', 'c_110100_b', 4, 'g--', 2);
figure(4);
subplot(221); legend('Theoretical', '111110 close', '111100 close', '011100 close', '110100 close b', '110100 close', '110100 close b');
subplot(223); legend('111110 close', '111100 close', '011100 close', '110100 close b', '110100 close', '110100 close b');
subplot(222); legend('Theoretical', '111110 close', '111100 close', '011100 close', '110100 close b', '110100 close', '110100 close b');
subplot(224); legend('111110 close', '111100 close', '011100 close', '110100 close b', '110100 close', '110100 close b');



% k_111111 = solve_disp(w, k, C_111111, '111111');
% k_111110 = solve_disp(w, k, C_111110, '111110');
% k_111100 = solve_disp(w, k, C_111100, '111100');
% k_111010 = solve_disp(w, k, C_111010, '111010');
% k_011100 = solve_disp(w, k, C_011100, '011100');
% k_110100 = solve_disp(w, k, C_110100, '110100');
% 
% k_c_111111 = solve_disp(w, k, C_c_111111, '111111');
% k_c_111110 = solve_disp(w, k, C_c_111110, '111110');
% k_c_111100 = solve_disp(w, k, C_c_111100, '111100');
% k_c_011100 = solve_disp(w, k, C_c_011100, '011100');
% k_c_011100_b = solve_disp(w, k, C_c_011100_b, '011100');
% k_c_110100 = solve_disp(w, k, C_c_110100, '110100');
% k_c_110100_b = solve_disp(w, k, C_c_110100_b, '110100');
