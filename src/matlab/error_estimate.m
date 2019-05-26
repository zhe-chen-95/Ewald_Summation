L = 2;
M = 53;
xi = 10;
m = 8;
P = 24;
pr = 2;
np = 1024;


C_F = 1; C_R = 50; 
eta = (P*L*xi / (M*m))^2;
E_F = C_F * exp(-pi^2 * M^2 / (4*L^2*xi^2));
E_Q = 4 * exp(-pi^2 * P^2 / (2*m^2)) + erfc(m / sqrt(2));
E_R = C_R * (1 / (xi^2) + pr / xi) * exp(-pr^2 * xi^2);

C_tR = 5e-7; C_tF = 7e-9; C_tG = 1e-7;
t_R = C_tR * np^2 * pr^3;
t_F = C_tF * (M^3) * log(M^3);
t_G = 2 * C_tG * np * P^3;
time = t_R + t_F + t_G + 0.5;
fprintf(['eta = ', num2str(eta), ' E_F = ', num2str(E_F), ' E_Q = ', num2str(E_Q), ' E_R = ', num2str(E_R), ' \n']);
fprintf(['time = ', num2str(time), ' t_F = ', num2str(t_F), ' t_G = ', num2str(t_G), ' t_R = ', num2str(t_R), ' \n']);