clear;clc

syndir = 'lab/';
% syn_list = {'GC_020.mat'; 'GC_025.mat'; 'GC_040.mat'; 'GC_050.mat'; 'GC_100.mat'};
% syn_list = {'GC_020.mat'; 'GC_025.mat'; 'GC_040.mat'; 'GC_050.mat'};
% syn_list = {'GC_020.mat'; 'GC_100_TZ.mat'; 'GC_100.mat'};
% syn_list = {'GC_025_TZ.mat'; 'GC_100_TZ.mat'};
% syn_list = {'GC_025_TZ.mat'; 'GC_100_TZ.mat'; 'GC_100.mat'; 'GC_050.mat'; 'GC_040.mat'; 'GC_025.mat'; 'GC_020.mat'};
% syn_list = {'GC_025_TZ.mat'; 'GC_100_TZ.mat'; 'GC_100.mat'; 'GC_050.mat'; 'GC_040.mat'; 'GC_025.mat'; 'GC_020.mat'};
syn_list = {'GC_020.mat', 'GC_025.mat', 'GC_100.mat', 'GC_025_TZ.mat', 'GC_100_TZ.mat'};

% n1 = 2;
% n2 = 2;

clim = 1;
i = 140;

clist = {'k', 'r', 'b', 'm', 'c', 'g', 'y'};

for l = 1 : length(syn_list)
    syn_file = [syndir syn_list{l}];
    dd = load(syn_file);
    if l == 1
        d0 = dd.d';
    end
    
    figure(1);
%     subplot(n1, n2, l);
    if l == 1
        imagesc(dd.d', [-clim clim]);
%     else
%         imagesc(dd.d' - d0, [-clim clim]);
    end
    
    figure(2);
    if l == 1
        plot(dd.t_axis, dd.d(i, :), [clist{l} '-'], 'linewidth', 2); hold on;
    else
        plot(dd.t_axis, dd.d(i, :), [clist{l} '--'], 'linewidth', 2); hold on;
    end
    
end

