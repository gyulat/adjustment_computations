% Sequential RANSAC line fitting
clear all; close all
% example data
d = load('star5.dat');
x = d(:,1);  y = d(:,2); 
tol = 0.05;  % fit threshold
k = 100;  % number of iterations
m = 5; % number of models
figure(1); hold on;
plot(x,y,"ko");
for n = 1:m
  nd = length(x);
  nmax = 0; % empty consensus set
  for i=1:k
    % select 2 random points
    is = randperm(nd,2);
    % determine parameters of line ax+by-1=0
    A = d(is,:); b = [1; 1]; p = A\b;
    % data distance from line
    t = abs(p(1)*x+p(2)*y-1)/sqrt(p(1)^2+p(2)^2);
    xk = x(t<tol);  yk = y(t<tol); % data that fit
    nin = length(xk);  % cardinality of consensus set
    if nin > nmax  % thus far the best
      xin = xk; yin = yk; nmax = nin;
      xout = x(t>=tol); yout = y(t>=tol); % not a consensus set
      bp = p;  % best fitting line
    end
  end
  pls = polyfit(xin,yin,1); % LSQ line fit
  plot(xin,yin,"go")
  line = @(x,y) bp(1)*x + bp(2)*y - 1;
  ezplot(line,[-1,1,-1,1]) % plot line
  % remove consensus set
  x = xout; y = yout; d = [x,y];
end
axis equal; ylim([-1,1]);
title("Sequential RANSAC line fitting");


