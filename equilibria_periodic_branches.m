clear;                           % clear variables
format compact
close all;                       % close figures

% % add path to ddebiftool (customise to your specific location)
% addpath('~/dde_biftool_v3.1.1/ddebiftool');  
% addpath('~/dde_biftool_v3.1.1/ddebiftool_utilities');  

addpath('/Users/ameliadodd/Documents/University Work/Year 4 /Mathematical Sciences Project/matlab/dde_biftool_v3.1.1/ddebiftool');  
addpath('/Users/ameliadodd/Documents/University Work/Year 4 /Mathematical Sciences Project/matlab/dde_biftool_v3.1.1/ddebiftool_utilities');  

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
    
% set up initial steady state 
stst.kind='stst';
stst.parameter=[1.1, 1, 0.9, 1.5, 0.5, 0.5, 1.5, 1.5, 0.5, 0];
stst.x=[4/15;7/15;4/15];

%% Create branch of non-trivial Equilibria

[eqbr,suc]=SetupStst(funcs,'contpar',ind_tau,'x',stst.x,'parameter',stst.parameter,...
    'max_bound',[ind_tau,1],'max_step',[ind_tau,0.01]);
if ~suc
    error('equilibrium not found');
end
figure(1);clf
eqbr=br_contn(funcs,eqbr,100);
xlabel('\tau');
ylabel('a_1')

%% Stability of equilibria
[eqnunst,dom,triv_defect,eqbr.point]=...
    GetStability(eqbr,'funcs',funcs,'points',2:length(eqbr.point));

%% Branch off at Hopf bifurcation

% stop script if there is no Hopf point
if isempty(find(diff(eqnunst)==2,1,'first')) == 1
      error('eqnunst: conditions for hopf bifurcation not met');
end

% initialise branch of periodic orbits
indhopf=find(diff(eqnunst)==2,1,'first');
[per,suc]=SetupPsol(funcs,eqbr,indhopf,'contpar',ind_tau,'degree',3,'intervals',50,...
    'print_residual_info',1);
if ~suc
    error('initialization of periodic orbits failed');
end

% continuation of periodic orbit in parameter
per.parameter.max_step=[0,0.02];
per=br_contn(funcs,per,100);
per.parameter.max_step=[];
per=br_contn(funcs,per,100);
per=br_stabl(funcs,per,0,1);

save('eq_per_branches.mat');
