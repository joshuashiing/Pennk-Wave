% 
% This is a forward modeling code to generate synthetic tested data as input
% for Q-RTM running using a viscoacoustic wave equation (Zhu and Harris, 2014). 
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
%% Article:
% ************************************************************
% Zhu T., and Harris J. M., (2014), 
% Modeling acoustic wave propagation in heterogeneous attenuating media using decoupled fractional Laplacians: 
% Geophysics, 79, no.3, T105-T116, doi:10.1190/geo2013-0245.1.
% ************************************************************
%


addpath ./misc

[sou,rec,para,vel,rho,qp] = para_FMmodel;


para.f_cutoff = 120;
para.taper = 0.2;

para.cq2sls = 1;
para.freq_ref_elastic = 1000;

tic
for izs = 1:size(sou,1)
    
    disp('shot #')
    izs
    
    % shot position [x z]
    sc=[sou(izs,1),sou(izs,2)]; 
    nt = para.nstep;
    para.sc=sc;para.izs = izs;
    
    % no data input
    syn_data = 0;
    
    % outupt data filename
    para.filename = ['Data/syn2Q',num2str(izs),'.mat'];
    % solving viscoacoustic wave equation
    FM_viscoacoustic(syn_data,vel,rho,qp,para)
    
end
toc;
%matlabpool close
