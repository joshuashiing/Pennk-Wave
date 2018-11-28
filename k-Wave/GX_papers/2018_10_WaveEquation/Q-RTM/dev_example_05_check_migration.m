clear;clc

migdir = 'Data_example_05/';
mig_list = {'mig_cd_rf.mat', 'mig_cd_ac.mat', 'mig_cd.mat', 'mig_cd_sq.mat'};
% mig_list = {'mig_cd_rf.mat', 'mig_cd_ac.mat'};

n1 = 4;
n2 = 1;
scale = 1;

cut = 1;
t_cut = 70;
% t_cut = 1;
clim = 5;

i = 220;
clist = {'k', 'r', 'b', 'm', 'c', 'g', 'y'};

% % imaging condition without laplacian
% for l = 1 : length(mig_list)
%     mig_file = [migdir mig_list{l}];
%     d = load(mig_file);
%     n = size(d.mig1, 3);
%     
%     img = sum(d.mig1(cut:end, :, :), 3);
%     img = del2(img);
%     
%     figure(1);
%     subplot(n1, n2, l);
%     clim = max(abs(img(:))) * scale;
%     imagesc(img, [-clim, clim]);
%     colorbar;
%     
%     figure(2);
%     plot(img(:, i), [clist{l} '-'], 'linewidth', 2); hold on;
% end

% imaging condition with laplacian
for l = 1 : length(mig_list)
    mig_file = [migdir mig_list{l}];
    d = load(mig_file);
    n = size(d.mig2, 3);
    
    img = sum(d.mig1(cut:end, :, :), 3);
    img = del2(img);
    
    figure(3);
    subplot(n1, n2, l);
%     clim = max(abs(img(:))) * scale;
    imagesc(img, [-clim, clim]);
    colorbar;
    
    figure(4);
    plot(img(:, i), [clist{l} '-'], 'linewidth', 2); hold on;
    xlim([t_cut, size(img, 1)]);
end
