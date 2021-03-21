function p_pplot_col(point,clr1,clr2,clr3)

% function p_pplot(point,component,colour)
% INPUT:
%       point point whose profile needs plotting
%	component optional component value
%	colour optional colour to plot with

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
if ~strcmp(point.kind,'psol') && ~strcmp(point.kind,'hcli')
  error('P_PPLOT: point type %s does not contain a profile.',point.kind);
end;

d=point.degree;
L=length(point.mesh);

if L==0
  L=size(point.profile,2);
  mesh=0:1/(L-1):1;
else
  mesh=point.mesh;
end;

hold on;
plot(mesh,point.profile(1,:),'color',clr1,'Linewidth',1.2);
plot(mesh,point.profile(2,:),'color',clr2,'Linewidth',1.2);
plot(mesh,point.profile(3,:),'color',clr3,'Linewidth',1.2);

return;


