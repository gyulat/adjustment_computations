% RANSAC line fitting minimal demo with plot
clc; clear all
% example data
d = load('linedata.txt');
x = d(:,1);  y = d(:,2); nd = length(x);
tol = 0.05;  % threshold for fit
k = 16;  % number of iterations
nmax = 0; % empty consensus set
figure(1)
for i=1:k
  % select two points at random
  is = randperm(nd,2);
  % plot of all/selected points
  plot(x,y,'marker','o','markerfacecolor','black')
  clf; hold on
  g = plot(x,y,'ko');
  set(g,'markerfacecolor','black');
  hold on
  h = plot(x(is),y(is),'ro');
  set(h,'markerfacecolor','red');
  % determine parameters of the line ax+by-1=0
  A = d(is,:); b = [1; 1]; p = A\b;
  % plot line
  line = @(x,y) p(1)*x + p(2)*y - 1;
  ezplot(line,[0,1,0,1]);
  ts = sprintf("iteration %d",i);
  title(ts);
  pause()
  % data distances to line
  t = abs(p(1)*x+p(2)*y-1)/sqrt(p(1)^2+p(2)^2);
  xk = x(t<tol);  yk = y(t<tol); % conform datum
  nin = length(xk);  % cardinality of the consensus set  
  if nin > nmax  % so far the best
    xin = xk; yin = yk; nmax = nin;
    bp = p;  % best line
  end
  % plot consensus set
  g = plot(xk,yk,'go');
  set(g,'markerfacecolor','green');
  ts = sprintf("iteration %d,  %d points fit",i,nin);
  title(ts);
  pause()
  hold off
end

% maximal consensus set
clf; hold on
g = plot(x,y,'ko');
set(g,'markerfacecolor','black');
h = plot(xin,yin,'go');
set(h,'markerfacecolor','green');
line = @(x,y) bp(1)*x + bp(2)*y - 1;
ezplot(line,[0,1,0,1]);
pls = polyfit(xin,yin,1); % LSQ line fit
bline = @(x) pls(1)*x + pls(2);
h = ezplot(bline,[0,1,0,1]);
set(h, 'Color','r')
ts = sprintf("cardinality of the consensus set: %d\n",nmax);
title(ts);
legend('outlier','conform','RANSAC', 'LSQ fit');
pause()
close all


