function FM_viscoacoustic(sensor_data,vel,rho,qp,para)
% To compare with acoustic analytic solution, 2D Constant Q modeling example
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

filename = para.filename;

%
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

% absorption properties
medium = medium_non_absorbing;
medium_CQ = medium_non_absorbing;

%%
Q = qp(1+PML_size:Ny+PML_size,1+PML_size:Nx+PML_size)';

% % beta in Carcione 2010's paper (Constant Q)
gama = 1/pi*atan(1./Q);

% old definition (Nov 19, 2012)
%medium_CQ.alpha_power = 1./(1.0 - gama); 

% New definition (Nov 19, 2012)
medium_CQ.alpha_model = 'CQ';
medium_CQ.alpha_power = gama;
medium_CQ.alpha_coeff = para.f0;      % Hz
medium_CQ.sound_speed = (para.f0/para.freq_ref_elastic).^gama.*medium.sound_speed;

%%
% store maximum supported frequency
f_max = kgrid.k_max*medium.sound_speed/(2*pi);

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
%t_max = para.t_max;
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
input_args = {'ReturnField',true,'Smooth', false, 'PMLInside', false, 'PMLAlpha',2,...
              'PMLSize', PML_size, 'PlotPML',false,'DataCast', 'single','PlotSim', false,...
              'DisplayMask','off', 'UsekSpace',false};

% run the forward simulation
 [sensor_z, wavefield_u]= RTM_SeisKspace_FW_TR_solver(kgrid, medium_CQ, source, sensor, input_args{:});
% 
 data_Q = sensor_z;
% % 
 save(filename,'data_Q','para')

end
