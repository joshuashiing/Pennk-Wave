clear;

% a = 1;
% b = -(6 + 1i);
% c = 11 + 5i;
% d = -(6 + 6i);

a = 1;
b = -6;
c = 11;
d = -6;


delta0 = b^2 - 3 * a * c;
delta1 = 2 * b^3 - 9 * a * b * c + 27 * a^2 * d;

cc = ((delta1 + (delta1^2 - 4 * delta0^3) ^ (1/2)) / 2) ^ (1/3);
xi1 = 1;
xi2 = -1 / 2 + sqrt(3) / 2 * 1i;
xi3 = -1 / 2 - sqrt(3) / 2 * 1i;
cc1 = cc * xi1;
cc2 = cc * xi2;
cc3 = cc * xi3;

x1 = -1 / (3*a) * (b + cc1 + delta0 / cc1);
x2 = -1 / (3*a) * (b + cc2 + delta0 / cc2);
x3 = -1 / (3*a) * (b + cc3 + delta0 / cc3);
