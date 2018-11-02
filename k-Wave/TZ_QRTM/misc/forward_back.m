function forward_back(izs,vel,rho,qp,syn_data,sc,rec,para)
%%
%%  Perform Forward / backward propagation
%%  put here due to parfor reason
%%  it is so slow to cross-product
%   so I will make it forward_back as a function to distribute to core
%   Tieyuan Zhu
%   APril 2012

filename = para.filename;
nt = para.nstep;
para.sc=sc;para.izs = izs;
%

mig = 0;mig_decon = 0;
epi = 1e-6;
s2 = 0;
s3 = 0;

%% with Q compensation
% Nov 19, 2012
[snapshot0,rtmsnapshot] = rtm_viscoacoustic(syn_data,vel,rho,qp,para);

nt = size(snapshot0,3);
tic
disp('mig xcorr')
for i = 1:nt-10
    
    % cross-correlation imaging condition
    mig = snapshot0(:,:,i).*rtmsnapshot(:,:,nt-i+1)+mig;
    s2 = snapshot0(:,:,i).^2+s2;
    s3 = rtmsnapshot(:,:,i).^2+s3;

end
toc

% save output parameters
save(filename, 'mig','s2','s3')

end