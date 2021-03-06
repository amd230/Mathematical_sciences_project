function br_plot_col(branch,x_m,y_m,clr)
%% plot branch
% function br_plot(branch,x_measure,y_measure,clr)
% INPUT:
%	branch branch of points
%	x_measure scalar measure for x-coordinate
%	y_measure scalar measure for y_coordinate
%	clr colour line to plot with

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
if ~exist('clr','var')
  clr='';
end;
if isempty(clr)
  clr='';
end;

ll=length(branch.point);

if ll<=1
  error('BR_PLOT: ll=%d, branch contains too few points.',ll);
end;

if ~isempty(x_m)
  [ffx,llx]=br_measr(branch,x_m);
else
  ffx=1:ll;
  llx=ones(ll);
end;
if ~isempty(y_m)
  [ffy,lly]=br_measr(branch,y_m);
else
  ffy=1:ll;
  lly=ones(ll);
end;

mx=max(llx);
my=max(lly);

if mx~=min(llx)
  for i=1:mx
    j_end=0;
    while 1
      llx(ll+1)=mx+1;
      j_start=j_end+1;
      while llx(j_start)<i
       j_start=j_start+1;
      end;
      if j_start==ll+1
        break;
      end;
      j_end=j_start+1;
      llx(ll+1)=0;
      while llx(j_end)>=i
        j_end=j_end+1;
      end;
      j_end=j_end-1;
      plot(ffx(j_start:j_end,i),ffy(j_start:j_end),'color',clr,'linewidth',1.2);
      hold on;
    end;
  end;
elseif my~=min(lly)
  for i=1:my
    j_end=0;
    while 1
      lly(ll+1)=my+1;
      j_start=j_end+1;
      while lly(j_start)<i
       j_start=j_start+1;
      end;
      if j_start==ll+1
        break;
      end;
      j_end=j_start+1;
      lly(ll+1)=0;
      while lly(j_end)>=i
        j_end=j_end+1;
      end;
      j_end=j_end-1;
      plot(ffx(j_start:j_end),ffy(j_start:j_end,i),'color',clr,'linewidth',1.2);
      hold on;
    end;
  end;
else
  plot(ffx,ffy,'color',clr,'linewidth',1.2); 
  hold on;
end;

return;
