clear;

% parameter values
s = [1.1,1,0.9];
ro = [0 1.5 0.5; 0.5 0 1.5; 1.5 0.5 0];
ddelay = 0.1; % vary this to change the value of tau
   
% define function for delayed system
ff = @(t,x,Z) [
    x(1)*(s(1)-x(1)-ro(1,2)*Z(2,1)-ro(1,3)*Z(3,1));
    x(2)*(s(2)-x(2)-ro(2,1)*Z(1,1)-ro(2,3)*x(3));
    x(3)*(s(3)-x(3)-ro(3,1)*Z(1,1)-ro(3,2)*x(2))];

% define function for non-delayed system
f = @(t,x,Z) [
    x(1)*(s(1)-x(1)-ro(1,2)*x(2)-ro(1,3)*x(3));
    x(2)*(s(2)-x(2)-ro(2,1)*x(1)-ro(2,3)*x(3));
    x(3)*(s(3)-x(3)-ro(3,1)*x(1)-ro(3,2)*x(2))];

% initial condition for non-delayed system
initial = [0.03,0.02,0.01];

% define end of time intervals for solutions
tt = 5000;
a = 10;

% solution of non-delayed system
sol = dde23(f,[],initial,[0,tt],ddeset('RelTol',1e-9,'AbsTol',1e-7));

% solution of non-delayed system to use as initial condition for the
% delayed case 
solh = dde23(f,[],initial,[0,a],ddeset('RelTol',1e-6));

% solution of delayed system
sol1 = dde23(ff,ddelay,solh,[a,tt],ddeset('RelTol',1e-9,'AbsTol',1e-9));

%% plot non-delayed solution
% (replace sol.x,sol.y with sol1.x,sol1.y to plot delayed solution - remember to change title)
figure(1); clf;
h=plot(sol.x,sol.y,'LineWidth',1.2);
set(h, {'color'}, {[.37 .65 .47];[.39 .58 .93];[1 .57 .69]}); % line colours
title('Trajectory Graph, $\tau = 0$','Fontsize',16,'Interpreter','latex');
legend('$a_1(t)$','$a_2(t)$','$a_3(t)$','Interpreter','latex','Fontsize',14);
xlabel('Time $t$','Interpreter','latex','Fontsize',14)
ylabel('$a_i(t)$','Interpreter','latex','Fontsize',14);
set(gca,'Fontsize',14,'FontName', 'CMU Serif'); % axis properties
ylim([0 1.2]);

%% plot delayed solution, with lines indicating equilibria values
figure(2); clf;
hold on;
hh=plot(sol1.x,sol1.y,'LineWidth',1.2);
h1=plot([0,tt],[1.1,1.1],':','LineWidth',1.2); % add equilibria markers
h2=plot([0,tt],[1,1],':','LineWidth',1.2);
h3=plot([0,tt],[0.9,0.9],':','LineWidth',1.2);
set(hh, {'color'}, {[.37 .65 .47];[.39 .58 .93];[1 .57 .69]});
set(h1, 'color', [.37 .65 .47]);
set(h2, 'color', [.39 .58 .93]);
set(h3, 'color', [1 .57 .69]);
hold off;
title(['Trajectory Graph, $\tau$ = ' num2str(ddelay)],'Fontsize',18,'Interpreter','latex');
legend('$a_1(t)$','$a_2(t)$','$a_3(t)$','Location','northwest','Interpreter','latex','Fontsize',16);
ylim([0 1.2]);
xlabel('Time $t$','Interpreter','latex','Fontsize',16)
ylabel('$a_i(t)$','Interpreter','latex','Fontsize',16);
set(gca,'Fontsize',16,'FontName', 'CMU Serif');

%% plot phase graph for delayed solution
% (replace sol1.x,sol1.y with sol.x,sol.y to plot non-delayed solution - remember to change title)
figure(3); clf
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
