function draw_ellipse(x,Sx,pconf,color,linewidth)
if (exist('linewidth')==0)
    linewidth=1; 
end;
if (exist('color')==0)
    color='black'; 
end;
t=0:0.01:2*pi; 
 p=x*ones(size(t))+sqrtm(-2*log(1-pconf)*Sx)*[cos(t);sin(t)];
 plot(p(1,:),p(2,:),color,'LineWidth',linewidth); 
end

 
