clear;clc

migfile = ['Data_example_01_HQ/mig.mat'];
% migfile = ['Data_example_01/mig.mat'];
load(migfile);

mig_s1 = zeros(size(mig1, 1), size(mig1, 2));
mig_s2 = zeros(size(mig2, 1), size(mig2, 2));

n = 30;

cut = 20;

figure(1);
for i = 1 : n
    mig_s1 = mig_s1 + mig1(:, :, i);
    img = mig1(cut:end, :, i);
    clim = max(abs(img(:)));
    imagesc(img, [-clim, clim]); colorbar;
    %pause
end
img = mig_s1(cut:end, :);
clim = max(abs(img(:)));
imagesc(img, [-clim, clim]); colorbar;
%pause

for i = 1 : n
    mig_s2 = mig_s2 + mig2(:, :, i);
    img = mig2(cut:end, :, i);
    clim = max(abs(img(:)));
    imagesc(img, [-clim, clim]); colorbar;
    %pause
end
img = mig_s2(cut:end, :);
clim = max(abs(img(:)));
imagesc(img, [-clim, clim]); colorbar;
%pause

%%
figure(3);
subplot(121);
img = mig_s1(cut:end, :);
clim = max(abs(img(:)));
imagesc(img, [-clim, clim]); colorbar;

subplot(122); 
% img = mig_s2(cut:end, :);
img = diff(mig_s1(cut:end, :),2,1);
clim = max(abs(img(:)));
imagesc(img, [-clim, clim]*0.1); colorbar;
