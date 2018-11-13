function plot_gather(t_axis, d, h, plim, figureid, subid)

offset_axis = (0 : (size(d, 1) - 1)) * h;
figure(figureid);
subplot(subid);
imagesc(offset_axis, t_axis, d', plim);