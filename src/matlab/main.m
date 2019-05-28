%Zhe Chen
%Script to run ewald summation and check convergence
%%
close all;
clc;
clear;
L = 1;
M = 8; %nx
xi = 1.2;
m = 8; %std of sigma of gaussian kernel
P = M; %points within support of Gaussian kernel
pr = 4; % layers of real space
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
M = 8; %nx
% xi = 1;
xi_all=[0.5,0.8,1];
m = 8; %std of sigma of gaussian kernel
P = M; %points within support of Gaussian kernel
% pr = 1; % layers of real space
np = 100; % # of points
numthreads = 4;
output="../../results/";
% pr_all=[2,4,8,16];
pr_all=linspace(2,10,5);
vel_rec=cell(length(pr_all),length(xi_all));

for i=1:length(pr_all)
    pr=pr_all(i);
    for j=1:length(xi_all)
        xi=xi_all(j);
        eta=(P*L/M*xi/m)^2;
        status = system(['./main ',num2str(numthreads),' ',num2str(M),' ',num2str(np),' ',num2str(P/2),' ',num2str(pr),' ',...
            num2str(L),' ',num2str(xi),' ',num2str(eta)],'-echo');
        if status~=0
            error("not successful to do ./main");
        end
        data=read_ewald('../../results/realspace.txt');
        vel_rec{i,j}=data(:,2:end);
    end
end
%%
%real space error
close all;
error=zeros(length(pr_all)-1,length(xi_all));
for j=1:length(xi_all)
    data0=vel_rec{end,j};
    for i=1:length(pr_all)-1
        data1=vel_rec{i,j};
        err=sqrt(sum((data1-data0).^2,2));
%         error(i,j)=norm(err,2)/sqrt(sum(sum(data0.^2,2)));
        error(i,j)=mean(err)/np/sqrt(3);
        %     data0=data1;
    end 
end
err_bound=@(pr,xi)(12*pi/xi^2+16*sqrt(pi)*pr/xi).*exp(-pr.^2.*xi^2);
figure;
hold on;
for j=1:length(xi_all)
    xi=xi_all(j);
    plot(pr_all(1:end-1),error,'color',[0 0.4470 0.7410]);
    plot(pr_all(1:end-1),err_bound(pr_all(1:end-1),xi),'--','color',[0.8500 0.3250 0.0980]);
end
ylabel('Relative Error of real space')
set(gca,'YScale','log','YLim',[1.0e-15,1.0e2],'FontSize',14)

%%
%k-space // nx
clc
cd ~/Cecil/Course/1st_year_spring/hpc/HPC/Ewald_Summation/src/conv
addpath ../matlab/
DIM=3;
L = 1;
% M = 8; %nx
M_all=linspace(24,64,5);
% M_all=20;
% xi = 1;
xi_all=[5,10,15]; 
m = 20; %std of sigma of gaussian kernel
P = 64; %points within support of Gaussian kernel
pr = 1; % layers of real space
np = 100; % # of points
numthreads = 4;
output="../../results/";
% pr_all=[2,4,8,16];
% pr_all=linspace(2,10,5);

vel_rec=cell(length(M_all),length(xi_all));

for i=1:length(M_all)
    M=M_all(i);
    for j=1:length(xi_all)
        xi=xi_all(j);
        P=M;
        eta=(P*L/M*xi/m)^2;
        status = system(['./main ',num2str(numthreads),' ',num2str(M),' ',num2str(np),' ',num2str(P/2),' ',num2str(pr),' ',...
            num2str(L),' ',num2str(xi),' ',num2str(eta)]);
        if status~=0
            error("not successful to do ./main");
        end
        data=read_ewald('../../results/kspace.txt');
        vel_rec{i,j}=data(:,2:end);
    end
end
%%
%k-space nx error
close all;
error=zeros(length(M_all)-1,length(xi_all));
for j=1:length(xi_all)
    data0=vel_rec{end,j};
    for i=1:length(M_all)-1
        data1=vel_rec{i,j};
        err=sqrt(sum((data1-data0).^2,2));
        error(i,j)=mean(err)/sqrt(3);
%         error(i,j)=norm(err,2)/sqrt(sum(sum(data0.^2,2)));
        %     data0=data1;
    end
end
err_bound=@(M,P,xi) 2*L^2/sqrt(pi) * (2*sqrt(pi)*M+3*xi*L).*exp(-M.^2*(pi/xi/L)^2)+4*exp(-pi^2*P.^2/(2*m^2*L^2))+erfc(m/sqrt(2));
figure;
hold on;
for j=1:length(xi_all)
    xi=xi_all(j);
    plot(M_all(1:end-1),error,'color',[0 0.4470 0.7410]);
    plot(M_all(1:end-1),err_bound(M_all(1:end-1),M_all(1:end-1),xi),'--','color',[0.8500 0.3250 0.0980]);
end
ylabel('Relative Truncation Error of k-space')
set(gca,'YScale','log','FontSize',14)