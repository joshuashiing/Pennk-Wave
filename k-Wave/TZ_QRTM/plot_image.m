%
% plot images from the results of RTM_main.m
% Tieyuan Zhu @ UT-Austin 2015
% @

clear all

addpath ./misc

%% Read in velocity model data and plot it

[sou,rec,para,vel,rho,qp] = para_FMmodel;
x = para.x;
z = para.z;
nx = length(x);
nz = length(z);
figure

load('rtm_images/rtmig1.mat');
% mig = mig;
Stacked1 = zeros(size(mig')); 
Stacked2 = zeros(size(mig'));
Stacked3 = zeros(size(mig')); 
n1=101;
n2=151;
n3=201;
mig = 0.0;
tic;
for is = 1:1:size(sou,1)
    
    %% read migration data
    load(['rtm_images/rtmig',num2str(is),'.mat']);
    mig = mig';s2 = s2'+1e0;
    Stacked1 = (Stacked1+mig);

end
% applying Laplacian filter to remove low-frequency noise
rtm_mig1 = rtm_lap(Stacked1(1:end,1:end)')/20;

figure;
imagesc(x,z(1:end),rtm_mig1(:,1:end)');axis image
xlabel('Distance (m)'); ylabel('Depth (m)');
colormap(flipud(gray));
make_axes_bold(gca);
