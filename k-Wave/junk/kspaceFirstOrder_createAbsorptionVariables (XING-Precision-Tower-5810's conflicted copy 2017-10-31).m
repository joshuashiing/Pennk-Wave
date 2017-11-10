% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to create
%     the absorption variables.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 26th November 2010
%     last update - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby

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

% define the lossy derivative operators and proportionality coefficients
if strcmp(equation_of_state, 'absorbing')
            
    % make sure the operators are positive and real
    medium.alpha_coeff = abs(real(medium.alpha_coeff));
    medium.alpha_power = abs(real(medium.alpha_power));
    
    % convert the absorption coefficient to nepers.(rad/s)^-y.m^-1
    medium.alpha_coeff = db2neper(medium.alpha_coeff, medium.alpha_power);

    % compute the absorbing fractional Laplacian operator and coefficient
    if ~(isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'no_absorption'))
        absorb_nabla1 = (kgrid.k).^(medium.alpha_power - 2); 
        absorb_nabla1(isinf(absorb_nabla1)) = 0;
        absorb_nabla1 = ifftshift(absorb_nabla1);
        absorb_tau = -2 .* medium.alpha_coeff .* medium.sound_speed.^(medium.alpha_power - 1);
    else
        absorb_nabla1 = 0;
        absorb_tau = 0;
    end
       
    % compute the dispersive fractional Laplacian operator and coefficient
    if ~(isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'no_dispersion'))
        absorb_nabla2 = (kgrid.k).^(medium.alpha_power-1); 
        absorb_nabla2(isinf(absorb_nabla2)) = 0;
        absorb_nabla2 = ifftshift(absorb_nabla2);            
        absorb_eta = 2 .* medium.alpha_coeff .* medium.sound_speed.^(medium.alpha_power) .* tan(pi .* medium.alpha_power / 2);
    else
        absorb_nabla2 = 0;
        absorb_eta = 0;
    end
        
    % pre-filter the absorption parameters if alpha_filter is defined (this
    % is used for time-reversal photoacoustic image reconstruction
    % with absorption compensation)
    if isfield(medium, 'alpha_filter')
                
        % update command line status
        disp('  filtering absorption variables...');        
        
        % frequency shift the absorption parameters
        absorb_nabla1 = fftshift(absorb_nabla1);
        absorb_nabla2 = fftshift(absorb_nabla2);
                        
        % apply the filter
        absorb_nabla1 = absorb_nabla1 .* medium.alpha_filter;
        absorb_nabla2 = absorb_nabla2 .* medium.alpha_filter;

        % shift the parameters back
        absorb_nabla1 = ifftshift(absorb_nabla1);
        absorb_nabla2 = ifftshift(absorb_nabla2); 
           
    end    

    % modify the sign of the absorption operators if alpha_sign is defined
    % (this is used for time-reversal photoacoustic image reconstruction
    % with absorption compensation)
    if isfield(medium, 'alpha_sign')
       if numel(medium.alpha_sign) == 2
           
           % if two parameters are given, apply the first to the absorption
           % parameter and the second to the disperion parameters
           absorb_tau = sign(medium.alpha_sign(1)) .* absorb_tau;
           absorb_eta = sign(medium.alpha_sign(2)) .* absorb_eta;
           
       else
           error('medium.alpha_sign must be given as a 2 element array controlling absorption and dispersion, respectively.');
       end
    end
% GXTEST    
elseif strcmp(equation_of_state, 'absorbing_TZ14')
    gamma = atan(1 ./ medium.Q) / pi;
    gamma_ref = mean(gamma(:));
    c0 = medium.sound_speed;
    c = c0 .* cos(pi * gamma / 2);
    w0 = medium.f0 * 2 * pi;
    absorb_eta = -c0.^(2*gamma) .* w0.^(-2*gamma) .* cos(pi*gamma);
    absorb_tau = -c0.^(2*gamma-1) .* w0.^(-2*gamma) .* sin(pi*gamma);
    absorb_nabla1 = (kgrid.k) .^ (gamma_ref*2 - 1);
    absorb_nabla1(isinf(absorb_nabla1)) = 0;
    absorb_nabla1 = ifftshift(absorb_nabla1);
    absorb_nabla2 = (kgrid.k) .^ (gamma_ref*2);
    absorb_nabla2(isinf(absorb_nabla2)) = 0;
    absorb_nabla2 = ifftshift(absorb_nabla2);
    clear c0 w0 gamma
elseif strcmp(equation_of_state, 'absorbing_TF17')
    % Calculate media parameter matrices
    gamma = atan(1 ./ medium.Q) / pi;
    c0 = medium.sound_speed;
    c = c0 .* cos(pi * gamma / 2);
    w0 = medium.f0 * 2 * pi;
    
    % unroll the media parameter matrices into vectors
    gamma_ur = gamma(:)';
    c_ur = c(:)';
    w0_ur = w0(:)';
    
    % Initialize B matrix column by column
    B_col1 = zeros(6, size(gamma_ur, 2));
    B_col2 = B_col1;
    B_col3 = B_col1;
    B_col4 = B_col1;
    B_col5 = B_col1;
    B_col6 = B_col1;
    
    % Ensemble B matrix column by column
    B_tmp = [ones(size(gamma_ur));
            (1 - gamma_ur); 
            (-1/2*gamma_ur + 1/2*gamma_ur.^2)];
    B_col1 = [B_tmp .* cos(-pi*gamma_ur/2);
              B_tmp .* sin(-pi*gamma_ur/2)];
    B_tmp = [ones(size(gamma_ur));
            (2 - 2*gamma_ur); 
            (1 - 3*gamma_ur + 2*gamma_ur.^2)];
    B_col2 = [B_tmp .* cos(-pi*gamma_ur);
              B_tmp .* sin(-pi*gamma_ur)];
    B_tmp = [ones(size(gamma_ur));
            (3 - 3*gamma_ur); 
            (3 - 15/2*gamma_ur + 9/2*gamma_ur.^2)];
    B_col3 = [B_tmp .* cos(-3/2*pi*gamma_ur);
              B_tmp .* sin(-3/2*pi*gamma_ur)];
    B_tmp = [ones(size(gamma_ur));
            (2 - gamma_ur); 
            (1 - 3/2*gamma_ur + 1/2*gamma_ur.^2)];
    B_col4 = [B_tmp .* cos(pi/2 - pi*gamma_ur/2);
              B_tmp .* sin(pi/2 - pi*gamma_ur/2)];
    B_tmp = [ones(size(gamma_ur));
            (3 - 2*gamma_ur); 
            (3 - 5*gamma_ur + 2*gamma_ur.^2)];
    B_col5 = [B_tmp .* cos(pi/2 - pi*gamma_ur);
              B_tmp .* sin(pi/2 - pi*gamma_ur)];
    B_tmp = [ones(size(gamma_ur));
            (4 - 3*gamma_ur); 
            (6 - 21/2*gamma_ur + 9/2*gamma_ur.^2)];
    B_col6 = [B_tmp .* cos(pi/2 - 3/2*pi*gamma_ur);
              B_tmp .* sin(pi/2 - 3/2*pi*gamma_ur)];
          
    % Solve B * A = [1; 2; 1; 0; 0; 0] and assign k coefficients
    C_k1_ur = zeros(size(gamma_ur));
    C_k2_ur = C_k1_ur;
    C_k3_ur = C_k1_ur;
    C_k4_ur = C_k1_ur;
    C_k5_ur = C_k1_ur;
    C_k6_ur = C_k1_ur;
    for i_grid = 1 : size(gamma_ur)
        B_grid = [B_col1(:, i_grid), B_col2(:, i_grid), ...
            B_col3(:, i_grid), B_col4(:, i_grid), ...
            B_col5(:, i_grid), B_col6(:, i_grid)];
        A_grid = B_grid \ [1; 2; 1; 0; 0; 0];
        C_k1_ur(i_grid) = A_grid(1) * c_ur(i_grid) * w0_ur(i_grid);
        C_k2_ur(i_grid) = A_grid(2) * c_ur(i_grid)^2;
        C_k3_ur(i_grid) = A_grid(3) * c_ur(i_grid)^3 / w0_ur(i_grid);
        C_k4_ur(i_grid) = A_grid(4) * c_ur(i_grid);
        C_k5_ur(i_grid) = A_grid(5) * c_ur(i_grid)^2 / w0_ur(i_grid);
        C_k6_ur(i_grid) = A_grid(6) * c_ur(i_grid)^3 / (w0_ur(i_grid)^2);
    end
    
    % Reshape the k coefficients
    absorb_C_k1 = reshape(C_k1_ur, size(gamma));
    absorb_C_k2 = reshape(C_k2_ur, size(gamma));
    absorb_C_k3 = reshape(C_k3_ur, size(gamma));
    absorb_C_k4 = reshape(C_k4_ur, size(gamma));
    absorb_C_k5 = reshape(C_k5_ur, size(gamma));
    absorb_C_k6 = reshape(C_k6_ur, size(gamma));
    
    clear *_ur B_col*
    
    absorb_nabla1 = (kgrid.k) .^ (-1/2);
    absorb_nabla1(isinf(absorb_nabla1)) = 0;
    absorb_nabla1 = ifftshift(absorb_nabla1);
    absorb_nabla2 = (kgrid.k) .^ (1/2);
    absorb_nabla2(isinf(absorb_nabla2)) = 0;
    absorb_nabla2 = ifftshift(absorb_nabla2);
    
% GXTEST    
    
end