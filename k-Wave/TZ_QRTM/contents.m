% Contents
%
% two main codes: forward modeling code and RTM imaging code
%% Forward modeling code to generate synthetic data, 
%  all common shot gather data will be saved by shot index #.
%  all paramters are controlled in para_FMmodel.m
%  
%  When you start, please open
%  Forwardmodeling_main.m 
%
%  run it and see if you have problems. Once you are good, please check
%  para_FMmodel.m to setup your model
%
%  when you run it, you may want to look at wave propagation, please set
%  parameter "PlotSim" to true in the subroutine FM_viscoacoustic.m line 122
%
%  para_FMmodel.m
%  FM_viscoacoustic.m

%% RTM code to migrate synthetic data to map the reflectors
% RTM_main.m 
% para_FMmodel.m
% forward_back.m
% RTM_filter.m
% rtm_viscoacoustic.m
%
% RTM_main.m is the main program to run.
% rtm_viscoacoustic.m is the core to complete forward and backward propagation
% forward_back.m is applying the image condition for calculated wavefields
%