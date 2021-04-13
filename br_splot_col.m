function br_splot_col(br,xm,ym,distance)

i_start=1;

while i_start<=length(br.point)
  s=p_dststb(br.point(i_start));
  i_end=i_start;
  while p_dststb(br.point(i_end))==s
    i_end=i_end+1;
    if i_end>length(br.point) 
        break; 
    end
  end
  if i_end>i_start+1
    b=br;
    b.point=b.point(i_start:i_end-1);
    if s==1, br_plot_col(b,xm,ym,[.31 .78 .47]); else br_plot_col(b,xm,ym,[.88 .24 .19]); end;
    hold on;
  end;
  i_start=i_end;
end;

for i=1:length(br.point)-1
  [s1,d1,c1]=p_dststb(br.point(i));
  [s2,d2,c2]=p_dststb(br.point(i+1));
  x1=p_measur(br.point(i),xm);
  x2=p_measur(br.point(i+1),xm);
  y1=p_measur(br.point(i),ym);
  y2=p_measur(br.point(i+1),ym);
  x=(abs(d1)*x2+abs(d2)*x1)/(abs(d1)+abs(d2));
  y=(abs(d1)*y2+abs(d2)*y1)/(abs(d1)+abs(d2));
  if s1~=s2
    if s1==1 
        plot([x1 x],[y1 y],'color',[.31 .78 .47],'linewidth',1.2); 
        plot([x x2],[y y2],'color',[.88 .24 .19],'linewidth',1.2);
    elseif s2==1
        plot([x1 x],[y1 y],'color',[.88 .24 .19],'linewidth',1.2); 
        plot([x x2],[y y2],'color',[.31 .78 .47],'linewidth',1.2); 
    else
        plot([x1 x2],[y1 y2],'color',[.88 .24 .19],'linewidth',1.2); 
    end
  end
  if s1~=s2
    if abs(d1)>abs(d2) 
        c1=c2; 
    end
    switch c1,
      case 0, plot(x,y,'o','MarkerSize',12,'linewidth',1.2,'col',[0.33 0.38 0.44]);
      case 1, plot(x,y,'x','MarkerSize',12,'linewidth',1.2,'col',[0.33 0.38 0.44]);
      case -1, plot(x,y,'s','MarkerSize',12,'linewidth',1.2,'col',[0.33 0.38 0.44]);
    end
  end
end
end
