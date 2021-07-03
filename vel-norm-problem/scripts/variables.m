% plot variable
n = 1000;
inc = (pi/4)/n;
phi = (0:1:n-1) * inc;

% constants
r_0 = 1;
a_T = 0.2;
mu = 1;

Beta = sqrt(4*mu/((3*pi+8)*a_T*r_0^2));

save('variables.mat');