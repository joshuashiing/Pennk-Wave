function [sou,rec,para,vel,rho,qp] = para_model

seis_survey = 1;

%% Read in velocity model data and plot it
load BP_cut_model_gas_chimney.mat
%load model_surface.mat

%smooth velp for RTM
% slowp = smoothn(1./velp,5);
% velp = 1./slowp;

velocityModel = velp;%figure;imagesc(velocityModel)

%% Q model
qp = Qp;
rho = 2200;

%%
[nz,nx] = size(velocityModel);

vmax = max(velocityModel(:));
vmin = min(velocityModel(:));

dx = 12.5;
dz = 12.5;
x = (1:nx)*dx-dx;
z = (1:nz)*dz-dz;

figure(11)
subplot(2,1,1)
imagesc(x*1e-3,z*1e-3,velocityModel*1e-3)
xlabel('Distance (km)'); ylabel('Depth (km)');
colorbar 'SouthOutside'
%title('Velocity Model');
axis image
set(gca,'XaxisLocation','top')
make_axes_bold(gca)

subplot(2,1,2)
imagesc(x*1e-3,z*1e-3,smoothn(qp,1));
colorbar 'SouthOutside'
xlabel('Distance (km)'); ylabel('Depth (km)');
axis image
set(gca,'XaxisLocation','top')
make_axes_bold(gca)
colormap parula

title('Q Model');
% hold on
% hshot = plot(x(1),z(1),'w*');
% hold off
% colormap(seismic)

%% Create shot gathers
% Use the velocity model to simulate a seismic survey.  The wave equations
% is solved using finite differences for a defined initial wavefield.

% calculate time step dt from stability crierion for finite difference
% solution of the wave equation.
dt = 0.9*min(min(dz./velocityModel/sqrt(2)));
dt = 1*1e-3;

% determine time samples nt from wave travelime to depth and back to
% surface
vmin = min(velocityModel(:));
nt = round(sqrt((dx*nx)^2 + (dz*nx)^2)*2/vmin/dt + 1);
nt = 3201;
t  = (0:nt-1).*dt;

dt = 1.0e-3;
% add region around model for applying absorbing boundary conditions (20
% nodes wide)
V = [repmat(velocityModel(:,1),1,20) velocityModel repmat(velocityModel(:,end),1,20)];
V = V';
vel_V = [repmat(V(:,1),1,20) V repmat(V(:,end),1,20)];
V = vel_V';

%V(end+1:end+20,:) = repmat(V(end,:),20,1);

%
r=0.001;npml=20;
freesurface='no';
[ddx,ddy,vel,rho,qp]=padgridpml(velocityModel,rho,qp,r,npml,dx,freesurface);

[nz0,nx0] = size(vel);
x0 = (1:nx0)*dx;
z0 = (1:nz0)*dz;
qp = gaussian_smoother(qp,x0,z0,dx*2);

subplot(2,1,2)
imagesc(x,z,qp)
xlabel('Distance (m)'); ylabel('Depth (m)');
set(gca,'XaxisLocation','top');colorbar
axis image
colormap((gray))

make_axes_bold(gca)
% figure;imagesc(x,z,qp)
% Define frequency parameter for ricker wavelet
f  = 20;
%vel = velocityModel;
%% Generate shots and save to file and video

%vidObj = VideoWriter('../videos/FaultModelShots.avi');
%open(vidObj);
data = zeros(size(nt,nx));
%figure(gcf)
rec_x=(22:1:nx+19)';
sou_x = (21:10:nx+20)';

hold on
plot((rec_x-20)*dx,dz,'k^',(sou_x-20)*dx,dz,'r*')

nr = length(rec_x);
nsou = length(sou_x);

rec = [rec_x, 25*ones(nr,1)];
sou = [sou_x, 25*ones(nsou,1)];

k=0;    
para.f0=f;para.dx=dx;para.dz=dz;para.dt=dt;
para.ddx=ddx;para.ddy=ddy;para.t_max = t(end);
para.nstep=para.t_max/dt;

%% Migration
% vidObj = VideoWriter('../videos/Onedippml_fd_RTM.avi');
% open(vidObj);

%colormap seismic %bone    
%sr=[(21:3:nz+20)',21*ones(nr,1)]; % shot position [x z]   
para.sr = [rec_x, 25*ones(nr,1)];
para.npml=20;
para.seis_survey=seis_survey;
para.x = x;para.z=z;
%% smooth velocity for migration
%vel=1./smooth2(1./vel,57);rho=1./smooth2(1./rho,57);qp=1./smooth2(1./qp,57);
%vel=1./smooth2(1./vel,5);rho=1./smooth2(1./rho,5);qp=1./smooth2(1./qp,5);