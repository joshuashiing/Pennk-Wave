clear;clc

migfile = 'Data_example_01/mig_cd_rf.mat';

scale = 0.1;

cut = 1;

d = load(migfile);
n = size(d.mig1, 3);

for i = 1 : n
    img = d.mig1(cut:end, :, i);
    img = del2(img);
    clim = max(abs(img(:))) * scale;
    imagesc(img, [-clim, clim]);
    colorbar;
    pause
end

figure(1);
img = sum(d.mig1(cut:end, :, :), 3);
clim = max(abs(img(:))) * scale;
imagesc(img, [-clim, clim]);
colorbar;

figure(2);
img = sum(d.mig1(cut:end, :, :), 3);
img = del2(img);
clim = max(abs(img(:))) * scale;
imagesc(img, [-clim, clim]);
colorbar;
