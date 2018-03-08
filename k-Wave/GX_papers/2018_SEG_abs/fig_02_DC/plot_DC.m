function plot_DC(casename, flag, figureid, plot_flag, lw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the pre-computed wavenumber for different cases and plot
% Author: Guangchi Xing
% Date: 03/03/2018

% clear;
% casename = 'case01';
% flag = '111111';
% figureid = 1;
% plot_flag = 'r:';

load([casename, '_k']);
w = 2 * pi * f;
cp = w ./ real(k);
alpha = -imag(k);

if ~strcmp(flag, '')
    k_sol = eval(['k_', flag]);
    cp_sol = w ./ real(k_sol);
    alpha_sol = -imag(k_sol);
    
    figure(figureid);
    
    subplot(221);
    plot(f, cp_sol, plot_flag, 'linewidth', lw);
    hold on;
    
    subplot(223);
    plot(f, (cp_sol - cp) ./ cp * 100, plot_flag, 'linewidth', lw);
    hold on;
    
    subplot(222);
    plot(f, alpha_sol, plot_flag, 'linewidth', lw);
    hold on;
    
    subplot(224);
    plot(f, (alpha_sol - alpha) ./ alpha * 100, plot_flag, 'linewidth', lw);
    hold on;
    
else
    figure(figureid);
    
    subplot(221);
    plot(f, cp, plot_flag, 'linewidth', lw);
    hold on;
    
    subplot(222);
    plot(f, alpha, plot_flag, 'linewidth', lw);
    hold on;
end


