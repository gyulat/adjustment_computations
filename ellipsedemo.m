% RANSAC ellipse fit demo (Octave 8.4.0/Matlab R2024a)
clear all; close all
% example data
d = load('ellipdata.txt'); % full ellipse
x = d(:,1);  y = d(:,2); nd = length(x);
tol = 0.25;  % fit threshold
k = 145;  % number of iterations
k = 2*k;  % increased by variance
nmax = 0; % empty consensus set
figure(1); hold on
plot(x,y,'ko')
for i=1:k
  while 1
    % select five random points
    is = randperm(nd,5);
    % x^2+Bxy+Cy^2+Dx+Ey+F=0 determine ellipse parameters
    A = [x(is).*y(is), y(is).^2, x(is), y(is), ones(5,1)];
    b = -x(is).^2; p = A\b;
    % ellipse, if B^2-4*C < 0  
    if (p(1)^2-4*p(2))<0
      break
    end
  end
  ell = @(x,y) x.^2+p(1)*x.*y+p(2)*y.^2+p(3)*x+p(4)*y+p(5);
  % gradient norm
  gnorm = @(x,y) sqrt((2*x+p(1)*y+p(3)).^2+(p(1)*x+2*p(2)*y+p(4)).^2);
  % gradient weighted fit residuals
  t = abs(ell(x,y))./gnorm(x,y);
  xk = x(t<tol);  yk = y(t<tol); % select fit
  nin = length(xk);  % cardinality of consensus set
  if nin > nmax  % so far the best
    xin = xk; yin = yk; nmax = nin;
    bp = p;  % best ellipse
  end
end

fprintf("cardinality of maximum consensus set: %d\n",nmax)
plot(xin,yin,'go')
ell = @(x,y) x.^2+bp(1)*x.*y+bp(2)*y.^2+bp(3)*x+bp(4)*y+bp(5);
ezplot(ell,[0,6,0,6])
% LSQ ellipse fit (Fitzgibbon et al. 1996)
D = [xin.*xin, xin.*yin, yin.*yin, xin, yin, ones(size(xin))];
S = D'*D;
C(6,6)=0; C(1,3)=2; C(2,2)=-1; C(3,1)=2;
[gevec, geval] = eig(inv(S)*C);
[posR, posC] = find(geval>0 & ~isinf(geval));
pls = gevec(:,posC);
a = (pls./pls(1)); a = a(2:end);
ells = @(x,y) x.^2+a(1)*x.*y+a(2)*y.^2+a(3)*x+a(4)*y+a(5);
h = ezplot(ells,[0,6,0,6]);
set(h,"Color","red");
title(["cardinality of consensus set: ",num2str(nmax)])
legend("data","conforming","ellipse","LSQ fit")
% center of ellipse a=(A B C D E F), Rosin(1999)
a = pls;
xc = (a(2)*a(5)-2*a(3)*a(4))/(4*a(1)*a(3)-a(2)^2);
yc = (a(2)*a(4)-2*a(1)*a(5))/(4*a(1)*a(3)-a(2)^2);
% semi major and minor axes
a1 = sqrt(-2*(a(6)-(a(3)*a(4)^2-a(2)*a(4)*a(5)+a(1)*a(5)^2)/(4*a(1)*a(3)-a(2)^2)) ...
            /(a(1)+a(3) - sqrt(a(2)^2+(a(1)-a(3))^2)) );
a2 = sqrt(-2*(a(6)-(a(3)*a(4)^2-a(2)*a(4)*a(5)+a(1)*a(5)^2)/(4*a(1)*a(3)-a(2)^2)) ...
            /(a(1)+a(3) + sqrt(a(2)^2+(a(1)-a(3))^2)) );
% orientation angle
theta = 0.5*atan(a(2)/(a(1)-a(3)));
fprintf("Ellipse parameters: \n");
fprintf("xc: %.3f  yc: %.3f\n",xc,yc);  % xc: 2.694  yc: 1.875, original: (3,2)
fprintf("a1: %.3f  a2: %.3f\n",a1,a2);  % a1: 1.129  a2: 1.590, original: 1.5, 2.5
fprintf("theta: %.1f deg\n", theta*180/pi); % theta: 36.4, original: 30


% LSQ fit without RANSAC
% LSQ ellipse fit (Fitzgibbon et al. 1996)
D = [x.*x, x.*y, y.*y, x, y, ones(size(x))];
S = D'*D;
C(6,6)=0; C(1,3)=2; C(2,2)=-1; C(3,1)=2;
[gevec, geval] = eig(inv(S)*C);
[posR, posC] = find(geval>0 & ~isinf(geval));
pls = gevec(:,posC);
a = (pls./pls(1)); a = a(2:end);
ells = @(x,y) x.^2+a(1)*x.*y+a(2)*y.^2+a(3)*x+a(4)*y+a(5);
figure(2); hold on
plot(x,y,'ko')
h = ezplot(ells,[0,6,0,6]);
set(h,"Color","red");
title("LSQ fit to all points")
legend("data","LSQ fit")
% center of ellipse a=(A B C D E F), Rosin(1999)
a = pls;
xc = (a(2)*a(5)-2*a(3)*a(4))/(4*a(1)*a(3)-a(2)^2);
yc = (a(2)*a(4)-2*a(1)*a(5))/(4*a(1)*a(3)-a(2)^2);
% semi major and minor axes
a1 = sqrt(-2*(a(6)-(a(3)*a(4)^2-a(2)*a(4)*a(5)+a(1)*a(5)^2)/(4*a(1)*a(3)-a(2)^2)) ...
            /(a(1)+a(3) - sqrt(a(2)^2+(a(1)-a(3))^2)) );
a2 = sqrt(-2*(a(6)-(a(3)*a(4)^2-a(2)*a(4)*a(5)+a(1)*a(5)^2)/(4*a(1)*a(3)-a(2)^2)) ...
            /(a(1)+a(3) + sqrt(a(2)^2+(a(1)-a(3))^2)) );
% orientation angle
theta = 0.5*atan(a(2)/(a(1)-a(3)));
fprintf("Ellipse parameters: \n");
fprintf("xc: %.3f  yc: %.3f\n",xc,yc);  
fprintf("a1: %.3f  a2: %.3f\n",a1,a2);  
fprintf("theta: %.1f deg\n", theta*180/pi);

pause()


