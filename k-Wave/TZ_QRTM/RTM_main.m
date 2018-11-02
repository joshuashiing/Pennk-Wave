% 
% This is a Q-RTM code to perform reverse-time migration using 
% a viscoacoustic wave equation (Zhu and Harris, 2014). The theory of Q-RTM
% was published in Geophysics (Zhu et al., 2014). Its principle is to
% compensate for attenuation effects for both source and receiver
% wavefields before applying RTM imaging condition. 
%
% Author: Tieyuan Zhu (tieyuanzhu@gmail.com)
% Start date: March 2012
% Last date: September 30, 2015
% 
%% This code can be used with permission. It will be restricted for 
%% academia research and won't allow to be commericialized.
%%
%% contact: tieyuanzhu@gmail.com
%%
% Reference paper:
% ************************************************************
% Zhu T., Harris J. M., and Biondi B., (2014), 
% Q-compensated reverse time migration: 
% Geophysics, 79, no.3, S77-S87, doi:10.1190/geo2013-0344.1.
% ************************************************************
%


addpath ./misc

[sou,rec,para,vel,rho,qp] = para_FMmodel;

%% low-pass filter parameters
para.f_cutoff = 100;
para.taper = .2;
para.freq_ref_elastic = 1000; % reference frequency at elastic limit

tic
for izs = 1:size(sou,1)
    
    % read seismic data in .mat
    load(['Data/syn2Q',num2str(izs),'.mat'],'data_Q')

    syn_data = data_Q;

    % source position (grid)
    sc=[sou(izs,1),sou(izs,2)]; % shot position [x z]
     
    nt = para.nstep;
    para.sc=sc;para.izs = izs;
    % outupt RTM data filename
    para.filename = ['rtm_images/rtmigQ_norotate_f40_2Qdata_xcorr',num2str(izs),'.mat'];
    
    % forward and backward propagation
    forward_back(izs,vel,rho,qp,syn_data,sc,rec,para);

end
toc;