% DESCRIPTION:
%       subscript to create the absorption variables
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 26th November 2010
%       last update - 7th February 2012
%
%       author      - Tieyuan Zhu
%       date        - 26th April 2012
%       Can't solve ringing artifacts in simulation in Heterogeneous media
%       last update - 11th Oct 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2012 Bradley Treeby and Ben Cox

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% Ny = kgrid.Ny;
% Nx = kgrid.Nx;
% [kx0,kz0] = derivative1_2D(kgrid.k,kgrid.dx,kgrid.dy);
% [ky, kx] = ndgrid( kz0,kx0);
% k = sqrt(kx.^2 + ky.^2);


%% Kjartansson constant Q
% define the lossy derivative operators and proportionality coefficients
if strcmp(equation_of_state, 'absorbing')
            
    % make sure the operators are positive and real
    medium.alpha_coeff = abs(real(medium.alpha_coeff));
    medium.alpha_power = abs(real(medium.alpha_power));

    % compute the absorbing fractional Laplacian operator and coefficient
    if ~(isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'no_absorption'))
        if strcmp(medium.alpha_model,'CQ')
        %%%----------------------------------------------------------------      
        %%% fractional laplacian ------ loss
        %%% I use fourier method to compute it in wavenumber domain
        absorb_nabla1 = (kgrid.k).^(2.0*mean(medium.alpha_power(:))-1.0); 
        absorb_nabla1(isinf(absorb_nabla1)) = 0;
        absorb_nabla1 = ifftshift(absorb_nabla1);
        absorb_tau = (c./(2.0*pi*medium.alpha_coeff)).^(2.0*medium.alpha_power).*c.^(-1.0).*sin(medium.alpha_power*pi);         
%       %%%----------------------------------------------------------------               

        elseif strcmp(medium.alpha_model,'KVQ')

        absorb_nabla1 = (kgrid.k).^(medium.alpha_power-1); 
        absorb_nabla1(isinf(absorb_nabla1)) = 0;
        absorb_nabla1 = ifftshift(absorb_nabla1);
        absorb_tau = -(medium.alpha_coeff.*c).^medium.alpha_power.*c.^(- 1)*sin(medium.alpha_power*pi/2);
        
        absorb_nabla11 = (kgrid.k).^(2*medium.alpha_power-1); 
        absorb_nabla11(isinf(absorb_nabla11)) = 0;
        absorb_nabla11 = ifftshift(absorb_nabla11);
        absorb_tau1 = -(medium.alpha_power-0.5).*(medium.alpha_coeff.*c).^(2*medium.alpha_power).*...
                       c.^(- 1)*sin(medium.alpha_power*pi/2)*cos(medium.alpha_power*pi/2);
            
        end
    else
        absorb_nabla1 = 0;
        absorb_tau = 0;
    end
       
    % compute the dispersive fractional Laplacian operator and coefficient
    if ~(isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'no_dispersion'))
        if strcmp(medium.alpha_model,'CQ')
      
        %%%----------------------------------------------------------------      
        %%% fractional laplacian ------ dispersion
        %%% I use fourier method to compute it in wavenumber domain
             absorb_nabla2 = (kgrid.k).^(2.0*mean(medium.alpha_power(:))) ; 
             absorb_nabla2(isinf(absorb_nabla2)) = 0;
             absorb_nabla2 = ifftshift(absorb_nabla2);
             absorb_eta = (c./(2.0*pi*medium.alpha_coeff)).^(2.0*medium.alpha_power).*cos(medium.alpha_power*pi);              
        %%%----------------------------------------------------------------      
        elseif strcmp(medium.alpha_model,'KVQ')

        absorb_nabla2 = (kgrid.k).^(medium.alpha_power); 
        absorb_nabla2(isinf(absorb_nabla2)) = 0;
        absorb_nabla2 = ifftshift(absorb_nabla2);
        absorb_eta = -(medium.alpha_coeff.*c).^medium.alpha_power.*cos(medium.alpha_power*pi/2);     

        absorb_nabla22 = (kgrid.k).^(2*medium.alpha_power); 
        absorb_nabla22(isinf(absorb_nabla22)) = 0;
        absorb_nabla22 = ifftshift(absorb_nabla22);
        absorb_eta2 = -0.5*medium.alpha_power*(medium.alpha_coeff.*c).^(2*medium.alpha_power).*(cos(medium.alpha_power*pi/2))^2;   
        
            
        end
    else
        absorb_nabla2 = 0;
        absorb_eta = 0;
    end
        absorb_nabla3 = 1;
    % pre-filter the absorption parameters if alpha_filter is defined (this
    % is used for time-reversal photoacoustic image reconstruction
    % with absorption compensation)
    if isfield(medium, 'alpha_filter');
                
        % update command line status
        disp('  filtering absorption variables...');        
        
        % frequency shift the absorption parameters
        absorb_nabla1 = fftshift(absorb_nabla1);
        absorb_nabla2 = fftshift(absorb_nabla2);
                        
        % apply the filter
        absorb_nabla1 = absorb_nabla1.*medium.alpha_filter;
        absorb_nabla2 = absorb_nabla2.*medium.alpha_filter;

        % Debug in two days and found this problems of time reverse of constant Q
        % Because I need to substract these "dispersion" effects from
        % simulation (Nov 18, 2012)
        absorb_nabla3 = absorb_nabla3.*medium.alpha_filter; 
        
        
        % shift the parameters back
        absorb_nabla1 = ifftshift(absorb_nabla1);
        absorb_nabla2 = ifftshift(absorb_nabla2); 
        absorb_nabla3 = ifftshift(absorb_nabla3);   
        
    end    
    
    absorb_eps=1;
    % modify the sign of the absorption operators if alpha_sign is defined
    % (this is used for time-reversal photoacoustic image reconstruction
    % with absorption compensation)
    if isfield(medium, 'alpha_sign')
       if numel(medium.alpha_sign) == 2
           % if two parameters are given, apply the first to the absorption
           % parameter and the second to the disperion parameters
           absorb_tau = sign(medium.alpha_sign(1))*absorb_tau;
           absorb_eta = sign(medium.alpha_sign(2))*absorb_eta;
       elseif numel(medium.alpha_sign) == 3 
           absorb_tau = sign(medium.alpha_sign(1))*absorb_tau;
           absorb_eta = sign(medium.alpha_sign(2))*absorb_eta;
           absorb_eps = sign(medium.alpha_sign(3))*absorb_eps;
       else
           error('medium.alpha_sign must be given as a 2 element array controlling absorption and dispersion, respectively.');
       end
    end

elseif strcmp(equation_of_state, 'absorbing-SLS')
            
    % make sure the operators are positive and real
    qp = medium.qualityfactor;
    Lrelax = medium.Lrelax;
    f0 = source.f0;
    sum_tetp=zeros(size(qp));
    
    for k = 1:Lrelax
       tau(:,:,k) = (sqrt(1. + qp.^2)-1)./(qp*f0(k)*pi*2);
       tep(:,:,k) = (sqrt(1. + qp.^2)+1)./(qp*f0(k)*pi*2); 
       sum_tetp(:,:) = sum_tetp(:,:) + tep(:,:,k)./tau(:,:,k);
    end
    
    MR = rho0.*c.^2*Lrelax./sum_tetp;
    %c = sqrt(MR./rho0);
    % pre-filter the absorption parameters if alpha_filter is defined (this
    % is used for time-reversal photoacoustic image reconstruction
    % with absorption compensation)
    absorb_nabla1=1;
    if isfield(medium, 'alpha_filter');
                
        % update command line status
        disp('  filtering absorption variables...');        
                 
        % apply the filter
        absorb_nabla1 = medium.alpha_filter;           
    end
    
    absorb_tau = 1;
    % modify the sign of the absorption operators if alpha_sign is defined
    % (this is used for time-reversal photoacoustic image reconstruction
    % with absorption compensation)
    if isfield(medium, 'alpha_sign')
       if numel(medium.alpha_sign) == 1
           % if two parameters are given, apply the first to the absorption
           % parameter and the second to the disperion parameters
           absorb_tau = medium.alpha_sign(1);
       else
           error('medium.alpha_sign must be given as a 2 element array controlling absorption and dispersion, respectively.');
       end
    end
end
