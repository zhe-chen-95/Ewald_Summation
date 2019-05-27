%Zhe Chen
%Script to run ewald summation and check convergence
%%
L = 1;
M = 5; %nx
xi = 0.5;
m = 8; %std of sigma of gaussian kernel
P = M; %points within support of Gaussian kernel
pr = 1; % layers of real space
np = 100; % # of points
fprintf('------------------------------------------------------------\n')
C_F = 1; C_R = 50; 
eta = (P*L*xi / (M*m))^2;
E_F = C_F * exp(-pi^2 * M^2 / (4*L^2*xi^2)); % k-space // nx (14)
E_Q = 4 * exp(-pi^2 * P^2 / (2*m^2)) + erfc(m / sqrt(2)); % quadrature error // P (25)
E_R = C_R * (1 / (xi^2) + pr / xi) * exp(-pr^2 * xi^2); % real space (16)

C_tR = 5e-7; C_tF = 7e-9; C_tG = 1e-7;
t_R = C_tR * np^2 * pr^3;
t_F = C_tF * (M^3) * log(M^3);
t_G = 2 * C_tG * np * P^3;
time = t_R + t_F + t_G + 0.5;
fprintf(['eta = ', num2str(eta), ' E_F = ', num2str(E_F), ' E_Q = ', num2str(E_Q), ' E_R = ', num2str(E_R), ' \n']);fprintf('------------------------------------------------------------\n')
C_F = 1; C_R = 50; 
eta = (P*L*xi / (M*m))^2;
E_F = C_F * exp(-pi^2 * M^2 / (4*L^2*xi^2)); % k-space // nx (14)
E_Q = 4 * exp(-pi^2 * P^2 / (2*m^2)) + erfc(m / sqrt(2)); % quadrature error // P (25)
E_R = C_R * (1 / (xi^2) + pr / xi) * exp(-pr^2 * xi^2); % real space (16)

C_tR = 5e-7; C_tF = 7e-9; C_tG = 1e-7;
t_R = C_tR * np^2 * pr^3;
t_F = C_tF * (M^3) * log(M^3);
t_G = 2 * C_tG * np * P^3;
time = t_R + t_F + t_G + 0.5;
fprintf(['eta = ', num2str(eta), ' E_F = ', num2str(E_F), ' E_Q = ', num2str(E_Q), ' E_R = ', num2str(E_R), ' \n']);

%%
%real space
clc
cd ~/Cecil/Course/1st_year_spring/hpc/HPC/Ewald_Summation/src/conv
addpath ../matlab/
DIM=3;
L = 1;
M = 5; %nx
xi = .2;
m = 8; %std of sigma of gaussian kernel
P = M; %points within support of Gaussian kernel
% pr = 1; % layers of real space
np = 100; % # of points
numthreads = 4;
output="../../results/";
pr_all=[1,2,4,8,16];
vel_rec=cell(length(pr_all),1);
eta=(P*L/M*xi/m)
for i=1:length(pr_all)
    pr=pr_all(i);
    status = system(['./main ',num2str(numthreads),' ',num2str(M),' ',num2str(np),' ',num2str(P/2),' ',num2str(pr)],'-echo');
    if status~=0
        error("not successful to do ./main");
    end
    data=read_ewald('../../results/realspace.txt');
    vel_rec{i}=data(:,2:end);
end
%%
%real space error
error=zeros(length(pr_all)-1,1);
data0=vel_rec{end};
for i=1:length(pr_all)-1
    data1=vel_rec{i};
    err=sqrt(sum((data1-data0).^2,2));
    error(i)=norm(err,2)/sqrt(sum(sum(data0.^2,2)));
%     data0=data1;
end
loglog(pr_all(1:end-1),error);
