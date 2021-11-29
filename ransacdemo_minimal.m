% RANSAC line fitting minimal demo
clc; clear all
% example data
d = load('linedata.txt');
x = d(:,1);  y = d(:,2); nd = length(x);
tol = 0.05;  % threshold for fit
k = 16;  % number of iterations
nmax = 0; % empty consensus set
for i=1:k
  % select two points at random
  is = randperm(nd,2);
  % determine parameters of the line ax+by-1=0
  A = d(is,:); b = [1; 1]; p = A\b;
  % data distances to line
  t = abs(p(1)*x+p(2)*y-1)/sqrt(p(1)^2+p(2)^2);
  xk = x(t<tol);  yk = y(t<tol); % conform datum
  nin = length(xk);  % cardinality of the consensus set
  if nin > nmax  % so far the best
    xin = xk; yin = yk; nmax = nin;
    bp = p;  % best line
  end
end
fprintf("cardinality of the consensus set: %d\n",nmax)
pls = polyfit(xin,yin,1); % LSQ line fit
fprintf("Equation of line of LSQ fit: y = %.3f*x %+.3f\n",pls(1),pls(2))



