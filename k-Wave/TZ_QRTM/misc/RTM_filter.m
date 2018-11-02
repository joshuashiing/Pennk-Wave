function medium_CQ = RTM_filter(sensor_data,kgrid,dt,medium,medium_CQ,para)
% get the average frequency spectrum of the simulated sensor data
num_signals = length(sensor_data(:,1));
dt0 = 1/dt;
[as_f, as] = spectrum(sensor_data(1, :), dt0);
for index = 2:num_signals
    [as_f, sp] = spectrum(sensor_data(index, :), dt0);
    as = as + sp;
end
as = as/num_signals;

% compute the relative power spectrum
ps = log10(as.^2);
offset = max(ps(:));
ps = ps - offset;

% get the frequency spectrum of a single measurement
as_sing = spectrum(sensor_data(1,:), 1/dt);
ps_sing = log10(as_sing.^2) - offset;

% scale the frequency variable
[f_sc scale prefix] = scaleSI(as_f(end));

% define the cutoff frequency for the filter
f_cutoff = para.f_cutoff;
taper_ratio = para.taper;
% =========================================================================
% IMAGE RECONSTRUCTION WITH ABSORPTION COMPENSATION
% =========================================================================

% create the filter to regularise the absorption parameters
%medium_CQ.alpha_filter = getAlphaFilter_old(kgrid, medium, f_cutoff,'Rectangular',taper_ratio);
medium_CQ.alpha_filter = getAlphaFilter(kgrid, medium, f_cutoff,taper_ratio);

% %plot the amplitude spectrum
% figure;
% plot(as_f*scale, ps_sing, 'r-');
% ylabel('Power Spectrum [dB]');
% xlabel(['Frequency [' prefix 'Hz]']);
% hold on;
% plot(as_f*scale, ps, 'k-');
% ylim = [-8 0];
% plot([f_cutoff*scale f_cutoff*scale], [ylim(1) ylim(2)], 'k--');   
% %plot([f_max*scale f_max*scale], [ylim(1) ylim(2)], 'k--');  
% %set(gca, 'XLim', [0 15], 'Ylim', ylim);
% 
% % plot the filter, reducing the number of plotting points by a factor of 10
% figure;
% mesh(medium_CQ.alpha_filter(1:10:end, 1:10:end), 'EdgeColor', 'black');
% axis tight; 
% view([-29, 46]);