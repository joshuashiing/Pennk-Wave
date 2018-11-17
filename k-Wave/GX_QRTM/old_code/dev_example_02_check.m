clear;clc

migfile1 = ['Data_example_02/mig_co.mat'];
migfile2 = ['Data_example_02/mig_co_ac.mat'];
% migfile = ['Data_example_01/mig.mat'];
d1 = load(migfile1);
d2 = load(migfile2);

mig_s1 = zeros(size(d1.mig, 1), size(d1.mig, 2));
mig_s2 = zeros(size(d2.mig, 1), size(d2.mig, 2));

n1 = size(d1.mig, 3);
n2 = size(d2.mig, 3);

cut = 1;
clim = 5;

figure(1);
for i = 1 : n1
    mig_s1 = mig_s1 + d1.mig(:, :, i);
    img = d1.mig(cut:end, :, i);
%     clim = max(abs(img(:)));
%     imagesc(img, [-clim, clim]); colorbar;
    imagesc(del2(img), [-clim, clim]); colorbar;
    pause
end
img = mig_s1(cut:end, :);
% clim = max(abs(img(:)));
% imagesc(img, [-clim, clim]); colorbar;
imagesc(del2(img), [-clim, clim]); colorbar;
pause

for i = 1 : n2
    mig_s2 = mig_s2 + d2.mig(:, :, i);
    img = d2.mig(cut:end, :, i);
%     clim = max(abs(img(:)));
%     imagesc(img, [-clim, clim]); colorbar;
    imagesc(del2(img), [-clim, clim]); colorbar;
    pause
end
img = mig_s2(cut:end, :);
% clim = max(abs(img(:)));
% imagesc(img, [-clim, clim]); colorbar;
imagesc(del2(img), [-clim, clim]); colorbar;
pause

%%
figure(2);
subplot(121);
img = mig_s1(cut:end, :);
% clim = max(abs(img(:)));
% imagesc(img, [-clim, clim]); colorbar;
imagesc(del2(img), [-clim, clim]); colorbar;

subplot(122);
img = mig_s2(cut:end, :);
% img = diff(mig_s1(cut:end, :),2,1);
% clim = max(abs(img(:)));
% imagesc(img, [-clim, clim]); colorbar;
imagesc(del2(img), [-clim, clim]); colorbar;