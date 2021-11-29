% sequential RANSAC circle fitting
clear all; close all
% example data
d = load('circle5.dat');
x = d(:,1);  y = d(:,2); 
tol = 0.02;  % fit threshold
k = 100;  % number of iterations
m = 5; % number of models
figure(1); hold on;
plot(x,y,"k*"); 
axis equal
for n = 1:m
  nd = length(x);
  nmax = 0; % no consensus set
  for i=1:k
    % select 3 random points
    is = randperm(nd,3);
    % parameter determination of the circle (x-x0)^2 + (y-y0)^2 + R^2 = 0 
    A = [d(is,:),ones(3,1)]; b = -[d(is,1).^2+d(is,2).^2]; p = A\b;
    xc = -0.5*p(1); yc = -0.5*p(2);
    R = sqrt((p(1)^2+p(2)^2)/4-p(3));
    % distances of points from circle
    t = abs(sqrt((x-xc).^2+(y-yc).^2)-R);
    xk = x(t<tol);  yk = y(t<tol); % conform data
    nin = length(xk);  % cardinality of consensus set
    if nin > nmax  % so far the best
      xin = xk; yin = yk; nmax = nin;
      xout = x(t>=tol); yout = y(t>=tol); % outliers
      bp = p;  % best parameters
    end
  end
  % LSQ circle fit
  A = [xin,yin,ones(length(xin),1)]; b = -[xin.^2+yin.^2]; p = A\b;
  xc = -0.5*p(1); yc = -0.5*p(2);
  R = sqrt((p(1)^2+p(2)^2)/4-p(3));
  plot(xin,yin,"g*")
  circle = @(x,y) sqrt((x-xc).^2+(y-yc).^2)-R;
  ezplot(circle,[-0.5,1,-0.5,1.5]) % plot circle
  % remove consensus set
  x = xout; y = yout; d = [x,y];
end
axis equal;
title("Sequential RANSAC circle fitting");
pause()


