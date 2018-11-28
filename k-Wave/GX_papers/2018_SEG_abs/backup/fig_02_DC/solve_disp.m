function [k] = solve_disp(w, k_ref, C, flag)
% Solve for wavenumber numerically given specific coefficients of wave
% equation
% Author: Guangchi Xing
% Date: 03/03/2018

fprintf(['Solving for k_', flag, '\n']);
k = zeros(size(w));
for i = 1 : length(w)
    w_i = w(i);
    lhs = w_i ^ 2;
    syms kk;
    switch flag
        case '111111'
            rhs = C(1)*kk + C(2)*kk^2 + C(3)*kk^3 + ...
                C(4)*(1i*w_i)*kk + C(5)*(1i*w_i)*kk^2 + C(6)*(1i*w_i)*kk^3;
        case '111110'
            rhs = C(1)*kk + C(2)*kk^2 + C(3)*kk^3 + ...
                C(4)*(1i*w_i)*kk + C(5)*(1i*w_i)*kk^2;
        case '111100'
            rhs = C(1)*kk + C(2)*kk^2 + C(3)*kk^3 + C(4)*(1i*w_i)*kk;
        case '111010'
            rhs = C(1)*kk + C(2)*kk^2 + C(3)*kk^3 + C(4)*(1i*w_i)*kk^2;
        case '011100'
            rhs = C(1)*kk^2 + C(2)*kk^3 + C(3)*(1i*w_i)*kk;
        case '110100'
            rhs = C(1)*kk + C(2)*kk^2 + C(3)*(1i*w_i)*kk;
    end
    eqn = (lhs == rhs);
    k_root = vpasolve(eqn, kk, k_ref(i));
    k_root = double(k_root);
    [~, i_root] = min(abs(k_root - k_ref(i)));
    k(i) = k_root(i_root);
end