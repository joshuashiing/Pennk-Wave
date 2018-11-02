function [Fwavefield, Rwavefield]= rtm_viscoacoustic(sensor_data,vel,rho,qp,para)
% To compare with viscoacoustic analytic solution, 2D Constant Q modeling example
% with fractional Laplacian. 
% 
% *** notice **
% careful with w0 and c0. 
%
% author: Tieyuan Zhu
% date: 3th May 2012
% last update: 4th May 2012
%  
% This function is part of the Seis-Wave Toolbox
% Copyright (C) 2012 Tieyuan Zhu

% This file is part of SEIS-Wave. SEIS-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% Seis-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with Seis-Wave. If not, see <http://www.gnu.org/licenses/>. 


addpath(genpath('./misc/'));

% npml = 20;freesur='no';
% [dxpml,dypml,vel,qp,rho]=padgridpml(vel,qp',rho,0.01,npml,0.5,freesur);
[nx,nz] = size(vel);

% =========================================================================
% FORWARD SIMULATION
% =========================================================================
simulation_sizex = nz;
simulation_sizey = nx;
% define the size of the simulation grid and the PML

PML_size = para.npml;                       % [grid points]
x = para.x;                          % [m]
y = para.z;                          % [m]

% reduce the number of pixels in Nx and Nz by the size of the PML so that
% using 'PMLInside' set to false will still give the correct simulation
% size 
Nx = simulation_sizex - 2*PML_size;  % [grid points]    
Ny = simulation_sizey - 2*PML_size;  % [grid points]
dx = para.dx;                           % [m]
dy = para.dz;                             % [m]

% create the computational grid
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium*ones(Nx,Ny)
medium_non_absorbing.sound_speed = vel(1+PML_size:Ny+PML_size,1+PML_size:Nx+PML_size)';	% [m/s]
medium_non_absorbing.density = rho(1+PML_size:Ny+PML_size,1+PML_size:Nx+PML_size)';	% [kg/m^3]

% create a duplicate of the propagation medium structure and append the
% absorption properties
medium = medium_non_absorbing;
medium_CQ = medium_non_absorbing;
medium_Q = medium_non_absorbing;
%%
medium_Q.qualityfactor = qp(1+PML_size:Ny+PML_size,1+PML_size:Nx+PML_size)';
Q = medium_Q.qualityfactor;

% % beta in Carcione 2010's paper (Constant Q)
gama = 1/pi*atan(1./Q);

% old definition (Nov 19, 2012)
%medium_CQ.alpha_power = 1./(1.0 - gama); 

% New definition (Nov 19, 2012)
medium_CQ.alpha_power = gama;

medium_CQ.alpha_coeff = para.f0;      % Hz

medium_CQ.sound_speed = (para.f0/para.freq_ref_elastic).^gama.*medium.sound_speed;

%
% %Kelvin-Voigt model
%  medium_Q.alpha_power = 1.;  
% % 
% % % Carcione(2007) P72 definition of Quality factor
% medium_Q.alpha_coeff = 1/(2*pi*250*Q);


medium_CQ.alpha_model = 'CQ';
%medium_Q.alpha_mode = 'no_dispersion';


% store maximum supported frequency
f_max = kgrid.k_max*medium.sound_speed/(2*pi);
% for k=1:medium.Lrelax
% source.f0(k) = 2^k;
% end
% alpha0 = pi/(medium.sound_speed*medium.qualityfactor); % neper/(Hz m)
% alpha_coeff=neper2db(alpha0);
source.f0= para.f0;


% create initial pressure distribution using makeDisc
disc_x_pos = para.sc(1)-PML_size;    % [grid points]
disc_y_pos = para.sc(2)-PML_size;    % [grid points]

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
rec_z = para.sr'-PML_size;
cart_sensor_mask = [(rec_z(1,:)-round(Nx/2))*dx; (rec_z(2,:)-round(Ny/2))*dy];
sensor.mask = cart_sensor_mask;

% create the time array used for the simulation, with t_max defined using
% Huygens' principle to avoid artifact trapping in the reconstruction
nt = para.nstep;
dt = para.dt;
t_max = (nt-1)*dt;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed, 0.1, t_max);
kgrid.t_array = 0:dt:t_max;
% create Ricker input 
f0= para.f0;
source_strength=1;
t0 = 1/f0;
rw=rickerwavelet(f0,dt,nt,t0,1);
%plot(t,rw), xlabel('Time'), ylabel('Amplitude')
source.p = source_strength.*rw;
source.p_mask = zeros(Nx,Ny);
source.p_mask(disc_x_pos,disc_y_pos)=1;
%source.p_mode= 'analytic';

% set the input options, switching off the smoothing (the input has already
% been smoothed), setting the PML to be outside the defined grid, casting
% to 'single' to speed up the example, and switching off visualisation
input_args = {'ReturnField',true,'Smooth', true, 'PMLInside', false, 'PMLAlpha',3,...
              'PMLSize', PML_size, 'PlotPML',false,'DataCast', 'single','PlotSim', false,...
               'UsekSpace',false};
           
 medium_CQ = RTM_filter(sensor_data,kgrid,dt,medium,medium_CQ,para);
% % % reverse the sign of the absorption proportionality coefficient
 medium_CQ.alpha_sign = [-1, 1, 1];        % [absorption, dispersion]; 
% % 

% run the forward simulation
[sensor_z, wavefield_u]= RTM_SeisKspace_FW_TR_solver(kgrid, medium_CQ, source, sensor, input_args{:});

% =========================================================================
% BACKPROPAGATION RECONSTRUCTION WITH ABSORPTION COMPENSATION
% =========================================================================
% create the computational grid
kgrid_recon = makeGrid(Nx, dx, Ny, dy);


% remove the initial pressure field from the source structure
source = rmfield(source, 'p');

% create a continuous binary sensor mask with the same radius as the
% Cartesian sensor mask used in the forward simulation
%pixel_radius = round(sensor_radius/kgrid_recon.dx);
%binary_sensor_mask = makeCircle(kgrid_recon.Nx, kgrid_recon.Ny, floor(kgrid_recon.Nx/2) + 1, floor(kgrid_recon.Ny/2) + 1, pixel_radius);

binary_sensor_mask = cart_sensor_mask;

% assign the sensor mask to the sensor structure
sensor.mask = binary_sensor_mask;

% attach the original time array
% add data for delay source
t_max0 = t_max - t0; kgrid_recon.t_array = kgrid.t_array;
%[kgrid_recon.t_array, dt] = makeTime(kgrid, medium.sound_speed, [], t_max);
%sensor_data0 = zeros(size(sensor_data,1),size(sensor_data,2)+round(t0/dt));
%sensor_data0(:,1:size(sensor_data,2)) = sensor_data;
% interpolate the simulated sensor data onto the continuous binary sensor
% mask to remove any gaps and assign to the time reversal field
%sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data, binary_sensor_mask, binary_sensor_mask, 'linear');
sensor.time_reversal_boundary_data = sensor_data;

% run the time-reversal reconstruction using the non absorbing medium
%[p0_recon, wavefield_ru] = SeiskspaceFirstOrder2D_rev(kgrid_recon, medium_non_absorbing, source, sensor, input_args{:});
%[p0_recon, wavefield_ru] = RTM_SeisKspace_FW_TR_solver(kgrid_recon, medium_non_absorbing, source, sensor, input_args{:});
%save(['wavefield_For_CQ/RTwavefield_noQ',num2str(para.izs),'.mat'],'wavefield_ru')
%load(['wavefield_For_CQ/RTwavefield_noQ',num2str(para.izs),'.mat']) 
% save p0_recon_model3.mat p0_recon wavefield_ru -v7.3
%load p0_recon_gauss.mat
% =========================================================================
% CHOOSING THE CUTOFF FREQUENCY
% =========================================================================
%  medium_CQ = RTM_filter(sensor_data,kgrid,dt,medium,medium_CQ,para);
% % %medium_CQ = RTM_filter(sensor_data,dt,medium,medium_CQ);
% % % reverse the sign of the absorption proportionality coefficient
medium_CQ.alpha_sign = [-1, 1, 1];        % [absorption, dispersion];

% run the time-reversal reconstruction
[~, wavefield_ru_Q] = RTM_SeisKspace_FW_TR_solver(kgrid_recon, medium_CQ, source, sensor, input_args{:});

Fwavefield = reshape(wavefield_u.wavefield_p,kgrid.Nx,kgrid.Ny,size(wavefield_u.wavefield_p,2));
Rwavefield = reshape(wavefield_ru_Q.wavefield_p,kgrid.Nx,kgrid.Ny,size(wavefield_ru_Q.wavefield_p,2));

end
