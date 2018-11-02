function [sou,rec,para,vel,rho,qp] = para_RTMmodel

seis_survey = 1;

%% Q model

rho = 2200;


%% Read in velocity model data and plot it
load vel.mat
velocityModel = vel';
velocityModel = ones(100,300);
 
[nz,nx] = size(velocityModel);

velocityModel(1:50,1:nx) = 1800;
% 
velocityModel(51:100,1:nx) = 2600;
velocityModel(1:end,1:end) = 3000;
 
qp = (velocityModel)/100;
% qp(find(velocityModel == 1600)) = 50;
% qp(find(velocityModel == 2000)) = 20;
% qp(find(velocityModel == 2400)) = 150;
qp(find(velocityModel == 2600)) = 100;
qp(find(velocityModel == 3000)) = 30;

rho(1:50,1:nx) = 1800;
% 
rho(51:100,1:nx) = 2200;

vmax = max(velocityModel(:));
vmin = min(velocityModel(:));

dx = 10;
dz = 10;
x = (1:nx)*dx - dx;
z = (1:nz)*dz - dz;

% velocityModel = 1./smoothn(1./velocityModel,10);
% qp = smoothn(qp,10);
%   qp = gaussian_smoother(qp,x,z,10);
%   velocityModel = gaussian_smoother(velocityModel,x,z,10);
 
figure
imagesc(x,z,velocityModel*1e-3)
xlabel('Distance (m)'); ylabel('Depth (m)');
%title('Velocity Model');
axis image;colorbar
set(gca,'XAxisLocation','top');colormap(flipud(gray))
set(gca,'XTick',[0 500 1000 1500 2000]);
set(gca,'YTick',[0 500 1000]);

figure
imagesc(x,z,qp)
xlabel('Distance (m)'); ylabel('Depth (m)');
%title('Q Model');
axis image;colorbar
set(gca,'XAxisLocation','top');colormap(flipud(gray))
set(gca,'XTick',[0 500 1000 1500 2000]);
set(gca,'YTick',[0 500 1000]);

hold on
hshot = plot(x(1),z(1),'w*');
hold off

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
nt = 1801;

t  = (0:nt-1).*dt;

% dt = 0.8e-3;
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

% [nz0,nx0] = size(vel);
% x = (1:nx0)*dx;
% z = (1:nz0)*dz;
% qp = gaussian_smoother(qp,x,z,20);

% figure;imagesc(x,z,qp)
% Define frequency parameter for ricker wavelet
f  = 25;
%vel = velocityModel;
%% Generate shots and save to file and video

%vidObj = VideoWriter('../videos/FaultModelShots.avi');
%open(vidObj);
data = zeros(size(nt,nx));
%figure(gcf)
rec_x=(22:1:nx+19)';
sou_x = (21:10:nx+20)';


nr = length(rec_x);
nsou = length(sou_x);

rec = [rec_x, 25*ones(nr,1)];
sou = [sou_x, 25*ones(nsou,1)];

hold on
plot((rec_x-npml)*dx,(rec(:,2)-npml)*dz,'k^',(sou_x-npml)*dx,(sou(:,2)-npml)*dz,'r*')

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