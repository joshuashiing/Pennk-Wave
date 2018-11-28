function [C] = solve_coef(B, C_norm, term_vec)
% Calculate the least square solution for the wave equation coefficient
% linear system
% Author: Guangchi Xing
% Date: 03/03/2018

B = B(:, term_vec);
C_norm = C_norm(term_vec);
A = (B' * B) \ (B' * [1; 2; 1; 0; 0; 0]);
C = A .* C_norm;