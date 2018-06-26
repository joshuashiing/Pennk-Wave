function [a_filter] = GX_FDRPfilter(alpha, Nmax, h)

% clear;clc
% 
% alpha = 1;
% Nmax = 5;
% h = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Riesz potential operator 

N = Nmax - 1;

d = 2;

C_alpha = 2^alpha * gamma(alpha / 2 + d / 2) / (pi^(d/2) * gamma(-alpha/2));
% C_alpha = gamma(1 - alpha / 2) / (2^alpha * pi^(alpha/2) * gamma(alpha/2));

%%%% Kernel function %%%%
f = @(x, y) (x.^2 + y.^2) .^ (-alpha / 2 - 1);
% I = zeros(N, N);

% parfor i = 0 : N
%     for j = 0 : N
%         xmin = (i - 1) * h;
%         xmax = (i + 1) * h;
%         ymin = (j - 1) * h;
%         ymax = (j + 1) * h;
%         I(i + 1, j + 1) = ...
%             integral2(f, xmin, xmax, ymin, ymax, 'RelTol', 1e-16, 'AbsTol', 0);
%     end
% end

I = zeros(N+1, N+1);
parfor i = 1:(N+1)
    for j = 1:(N+1)
        % -----------------------------------------------------------------
        xmin = (i-2)*h; xmax = i*h; ymin = (j-2)*h; ymax = j*h;
        I(i,j) = integral2(f, xmin, xmax, ymin, ymax, 'RelTol', 1e-16, 'AbsTol', 0);
        % -----------------------------------------------------------------
    end
end


a = I * C_alpha / 4;
a_q1 = flipud(a(2:end, :));
a_q2 = fliplr(a_q1(:, 2:end));
a_q3 = fliplr(a(:, 2:end));
a_q4 = a;
a_filter = -[a_q2, a_q1; a_q3, a_q4];
