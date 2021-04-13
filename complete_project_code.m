clear;                           % clear variables
format compact
close all;                       % close figures

% % add path to ddebiftool (customise to your specific location)
% addpath('~/dde_biftool_v3.1.1/ddebiftool');  
% addpath('~/dde_biftool_v3.1.1/ddebiftool_utilities');  

addpath('/Users/ameliadodd/Documents/University Work/Year 4 /Mathematical Sciences Project/matlab/dde_biftool_v3.1.1/ddebiftool');  
addpath('/Users/ameliadodd/Documents/University Work/Year 4 /Mathematical Sciences Project/matlab/dde_biftool_v3.1.1/ddebiftool_utilities');  

%% Find equilibria of system by solving RHS equal to zero 

% parameter values
s = [1.1,1,0.9];
ro = [0 1.5 0.5; 0.5 0 1.5; 1.5 0.5 0];

syms x y z
eqn1 = x*(s(1)-x-ro(1,2)*y-ro(1,3)*z) == 0;
eqn2 = y*(s(2)-y-ro(2,1)*x-ro(2,3)*z) == 0;
eqn3 = z*(s(3)-z-ro(3,1)*x-ro(3,2)*y) == 0;
sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
mysol = [sol.x(:),sol.y(:),sol.z(:)]; % matrix of equilibria
 
% find Jacobian of system for each positive equilibria
syms x y z
J=jacobian([x*(s(1)-x-ro(1,2)*y-ro(1,3)*z), y*(s(2)-y-ro(2,1)*x-ro(2,3)*z), z*(s(3)-z-ro(3,1)*x-ro(3,2)*y)], [x, y, z]);
J1 = subs(J, {x,y,z}, {0,0,0});
J2 = subs(J, {x,y,z}, {1.1,0,0});
J3 = subs(J, {x,y,z}, {0,1,0});
J4 = subs(J, {x,y,z}, {0,0,0.9});
J5 = subs(J, {x,y,z}, {4/15,7/15,4/15});

% find eigenvalues and eigenvectors of equilibria 
[V,D] = eig(J1);
[V1,D1] = eig(J2);
[V2,D2] = eig(J3);
[V3,D3] = eig(J4);
[V4,D4] = eig(J5);

%% Simulation for tau = 0 

% define function for non-delayed system
f = @(t,x,Z) [
    x(1)*(s(1)-x(1)-ro(1,2)*x(2)-ro(1,3)*x(3));
    x(2)*(s(2)-x(2)-ro(2,1)*x(1)-ro(2,3)*x(3));
    x(3)*(s(3)-x(3)-ro(3,1)*x(1)-ro(3,2)*x(2))];

% initial condition for non-delayed system
initial = [0.03,0.02,0.01];

% define end of time interval for solution
tt = 5000;

% solution of non-delayed system
sol = dde23(f,[],initial,[0,tt],ddeset('RelTol',1e-9,'AbsTol',1e-7));

%% plot non-delayed solution
% trajectory plot
figure(1); clf;
h=plot(sol.x,sol.y,'LineWidth',1.2);
set(h, {'color'}, {[.37 .65 .47];[.39 .58 .93];[1 .57 .69]}); % line colours
title('Trajectory Graph, $\tau = 0$','Fontsize',16,'Interpreter','latex');
legend('$a_1(t)$','$a_2(t)$','$a_3(t)$','Interpreter','latex','Fontsize',14);
xlabel('Time $t$','Interpreter','latex','Fontsize',14)
ylabel('$a_i(t)$','Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',14,'FontName', 'CMU Serif'); % axis properties
ylim([0 1.2]);

% plot phase graph for delayed solution
figure(2); clf
hold all;
plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:), 'color',[.6 .4 .8],'LineWidth',1.2);
h=plot3(0,0,0,'ko'); % add equilibria markers
plot3(0,0,9/10,'ko')
plot3(0,1,0,'ko')
plot3(11/10,0,0,'ko')
plot3(4/15,7/15,4/15,'ko')
grid on
title('Phase Graph, $\tau$ = 0','Fontsize',16,'Interpreter','latex');
xlabel('$a_1(t)$','Interpreter','latex','Fontsize',16)
ylabel('$a_2(t)$','Interpreter','latex','Fontsize',16);
zlabel('$a_3(t)$','Interpreter','latex','Fontsize',16)
legend(h,'Equilibria','Fontsize',14,'Interpreter','latex');
xlim([0 1.2]);
ylim([0 1.2]);
zlim([0 1.2]);
xticks(0:0.2:1.2);
yticks(0:0.2:1.2);
zticks(0:0.2:1.2);
set(gca,'Fontsize',14,'FontName', 'CMU Serif');
hold off;
view(3)

%% Define funcs for DDE-BIFTOOL

% define location of parameters
% param = [sigma1, sigma2, sigma3, ro12, ro13, ro21, ro23, ro31, ro32, tau]
ind_sigma1=1; ind_sigma2=2; ind_sigma3=3; 
ind_ro12=4; ind_ro13=5; ind_ro21=6; ind_ro23=7; ind_ro31=8; ind_ro32=9;
ind_tau=10;

% set up RHS system function 
% xx = [x1(t), x1(t-tau); x2(t), x2(t-tau); x3(t) x3(t-tau)]
my_sys_rhs=@(xx,par)[
    xx(1,1)*(par(1)-xx(1,1)-par(4)*xx(2,2)-par(5)*xx(3,2));
    xx(2,1)*(par(2)-xx(2,1)-par(6)*xx(1,2)-par(7)*xx(3,1));
    xx(3,1)*(par(3)-xx(3,1)-par(8)*xx(1,2)-par(9)*xx(2,1))];
% define location of delay parameter
my_sys_tau=@()(ind_tau); 

% store user-defined functions about system in a structure
funcs=set_funcs(...
    'sys_rhs',my_sys_rhs,...
    'sys_tau',my_sys_tau);

% Import branches for equilibrium and periodic orbits 
load('eq_per_branches.mat');

%% Bifurcation diagram
figure(3);clf
[xm,ym]=df_measr(0,eqbr);
br_splot_col(eqbr,xm,ym);
title('Bifurcation Diagram - varying $\tau$','Interpreter','latex','Fontsize',18);
xlabel('$\tau$','Interpreter','latex','Fontsize',16)
ylabel('$a_1(t)$','Interpreter','latex','Fontsize',16)
[xm,ym]=df_measr(0,per,'psol');
ym.col='max';
br_splot_col(per,xm,ym);
ym.col='min';
br_splot_col(per,xm,ym);
[xm,ym]=df_measr(0,per1,'psol');
ym.col='max';
br_splot_col(per1,xm,ym);
ym.col='min';
br_splot_col(per1,xm,ym);
[xm,ym]=df_measr(0,per2,'psol');
ym.col='max';
br_splot_col(per2,xm,ym);
ym.col='min';
br_splot_col(per2,xm,ym);
set(gca,'Fontsize',16,'FontName', 'CMU Serif');

% zoom in to bifurcation diagram
figure(4);clf
[xm,ym]=df_measr(0,eqbr);
br_splot_col(eqbr,xm,ym);
xlim([0,1])
title('Bifurcation Diagram - varying $\tau$','Interpreter','latex','Fontsize',18);
xlabel('$\tau$','Interpreter','latex','Fontsize',16)
ylabel('$a_1(t)$','Interpreter','latex','Fontsize',16)
[xm,ym]=df_measr(0,per,'psol');
ym.col='max';
br_splot_col(per,xm,ym);
ym.col='min';
br_splot_col(per,xm,ym);
set(gca,'Fontsize',16,'FontName', 'CMU Serif');

%% Periodic along the branch of periodic orbits 
figure(5); clf;
[xm,ym]=df_measr(0,per);
ym.field='period';
ym.col=1;
bb = [.6 .4 .8];
br_plot_col(per,xm,ym,bb); % look at the period along the branch:
xlabel('$\tau$','Interpreter','latex','Fontsize',14);ylabel('Period','Interpreter','latex','Fontsize',14);
title('Period of Periodic Orbits Against $\tau$','Interpreter','latex','Fontsize',16);
set(gca,'Fontsize',14,'FontName', 'CMU Serif');

%% Simulations for tau = 0.06

% define function for delayed system
ff = @(t,x,Z) [
    x(1)*(s(1)-x(1)-ro(1,2)*Z(2,1)-ro(1,3)*Z(3,1));
    x(2)*(s(2)-x(2)-ro(2,1)*Z(1,1)-ro(2,3)*x(3));
    x(3)*(s(3)-x(3)-ro(3,1)*Z(1,1)-ro(3,2)*x(2))];

% define the delay 
ddelay = 0.06;

% solution of non-delayed system to use as initial condition for the
% delayed case 
solh = dde23(f,[],initial,[0,10],ddeset('RelTol',1e-6));

% solution of delayed system
sol1 = dde23(ff,ddelay,solh,[10,tt],ddeset('RelTol',1e-9,'AbsTol',1e-9));

% phase plot 
figure(6); clf
hold all;
plot3(sol1.y(1,:),sol1.y(2,:),sol1.y(3,:), 'color',[.6 .4 .8],'LineWidth',1.2);
plot3(0,0,0,'ko'); % add equilibria markers
plot3(0,0,9/10,'ko')
plot3(0,1,0,'ko')
plot3(11/10,0,0,'ko')
plot3(4/15,7/15,4/15,'ko')
grid on
title(['Phase Graph, $\tau$ = ' num2str(ddelay)],'Fontsize',16,'Interpreter','latex');
xlabel('$a_1(t)$','Interpreter','latex','Fontsize',16)
ylabel('$a_2(t)$','Interpreter','latex','Fontsize',16);
zlabel('$a_3(t)$','Interpreter','latex','Fontsize',16)
legend(h,'Equilibria','Fontsize',14,'Interpreter','latex');
xlim([0 1.2]);
ylim([0 1.2]);
zlim([0 1.2]);
xticks(0:0.2:1.2);
yticks(0:0.2:1.2);
zticks(0:0.2:1.2);
set(gca,'Fontsize',14,'FontName', 'CMU Serif');
hold off;
view(3)

%% simulations for tau = 0.67

% redefine the delay 
ddelay1 = 0.67;

% solution of delayed system
sol2 = dde23(ff,ddelay1,solh,[10,tt],ddeset('RelTol',1e-9,'AbsTol',1e-9));

% plot delayed solution, with lines indicating carrying capacities 
figure(7); clf;
hold on;
hh=plot(sol2.x,sol2.y,'LineWidth',1.2);
h1=plot([0,tt],[1.1,1.1],':','LineWidth',1.2); % add equilibria markers
h2=plot([0,tt],[1,1],':','LineWidth',1.2);
h3=plot([0,tt],[0.9,0.9],':','LineWidth',1.2);
set(hh, {'color'}, {[.37 .65 .47];[.39 .58 .93];[1 .57 .69]});
set(h1, 'color', [.37 .65 .47]);
set(h2, 'color', [.39 .58 .93]);
set(h3, 'color', [1 .57 .69]);
hold off;
title(['Trajectory Graph, $\tau$ = ' num2str(ddelay1)],'Fontsize',18,'Interpreter','latex');
legend('$a_1(t)$','$a_2(t)$','$a_3(t)$','Location','northwest','Interpreter','latex','Fontsize',16);
ylim([0 1.2]);
xlabel('Time $t$','Interpreter','latex','Fontsize',16)
ylabel('$a_i(t)$','Interpreter','latex','Fontsize',16);
set(gca,'Fontsize',16,'FontName', 'CMU Serif');

% plot phase graph for delayed solution
figure(8); clf
hold all;
plot3(sol2.y(1,:),sol2.y(2,:),sol2.y(3,:), 'color',[.6 .4 .8],'LineWidth',1.2);
plot3(0,0,0,'ko'); % add equilibria markers
plot3(0,0,9/10,'ko')
plot3(0,1,0,'ko')
plot3(11/10,0,0,'ko')
plot3(4/15,7/15,4/15,'ko')
grid on
title(['Phase Graph, $\tau$ = ' num2str(ddelay1)],'Fontsize',16,'Interpreter','latex');
xlabel('$a_1(t)$','Interpreter','latex','Fontsize',16)
ylabel('$a_2(t)$','Interpreter','latex','Fontsize',16);
zlabel('$a_3(t)$','Interpreter','latex','Fontsize',16)
legend(h,'Equilibria','Fontsize',14,'Interpreter','latex');
xlim([0 1.2]);
ylim([0 1.2]);
zlim([0 1.2]);
xticks(0:0.2:1.2);
yticks(0:0.2:1.2);
zticks(0:0.2:1.2);
set(gca,'Fontsize',14,'FontName', 'CMU Serif');
hold off;
view(3)

%% continuation of hopf point in tau and sigma 

[hbranch,suc]=SetupHopf(funcs,eqbr,indhopf,'contpar',[ind_tau,ind_sigma1],'dir',ind_tau,'step',1e-3);
figure(9);clf;
hbranch=br_contn(funcs,hbranch,60);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,60);

%
[xm,ym]=df_measr(0,hbranch);
figure(9);clf;
bb = [.6 .4 .8];
br_plot_col(hbranch,xm,ym,bb);
xlabel('$\tau$','Interpreter','latex','Fontsize',12);
ylabel('$\sigma_1$','Interpreter','latex','Fontsize',12);
title('Hopf Point Continued in $\tau$ and $\sigma_{1}$','Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12,'FontName', 'CMU Serif');

% compute lyapnov coefficients to find if Hopf bifurcation is supercritical (L1<0) or subcritical (L1>0)
[L1,L1low]=HopfLyapunovCoefficients(funcs,hbranch);
fprintf('maximal L1 coefficient along Hopf branch: %g\n',max(L1));
fprintf('max of error estimate for L1 coefficient: %g\n',norm(L1-L1low,'inf'));

%% continuation of hopf point in tau and rho_31

[hbranch1,suc]=SetupHopf(funcs,eqbr,indhopf,'contpar',[ind_tau,ind_ro31],'dir',ind_tau,'step',1e-4);

figure(11);clf
hbranch1=br_contn(funcs,hbranch1,250);
hbranch1=br_rvers(hbranch1);
hbranch1=br_contn(funcs,hbranch1,60);
[xm,ym]=df_measr(0,hbranch1);
%%
figure(10);clf;
bb = [.6 .4 .8];
br_plot_col(hbranch1,xm,ym,bb);
xlabel('$\tau$','Interpreter','latex','Fontsize',12);
ylabel('$\rho_{31}$','Interpreter','latex','Fontsize',12);
title('Hopf Point Continued in $\tau$ and $\rho_{31}$','Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',12,'FontName', 'CMU Serif');

% compute lyapnov coefficients to find if Hopf bifurcation is supercritical (L1<0) or subcritical (L1>0)
[L11,L11low]=HopfLyapunovCoefficients(funcs,hbranch1);
fprintf('maximal L1 coefficient along Hopf branch: %g\n',max(L11));
fprintf('max of error estimate for L1 coefficient: %g\n',norm(L11-L11low,'inf'));
