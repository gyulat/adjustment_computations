% Solve overdetermined system of 3 equations with Kalman filter
%
% 2*x1 + 3*x2 = 18 + r1
% 4*x1 + 3*x2 = 24 + r2
%   x1 - x2 = 0 + r2
%
% variances of r1,r2,r3 are 1,4,4, respectively

clc; clear all; close all

% plot the system with zero residuals
eq1 = @(x1,x2) 2*x1+3*x2-18;
eq2 = @(x1,x2) 4*x1+3*x2-24;
eq3 = @(x1,x2) x1-x2;
figure(1)
hold on
ezplot(eq1,[0,10])
ezplot(eq2,[0,10])
ezplot(eq3,[0,10])
legend('2*x1+3*x2=18','4*x1+3*x2=24','x1-x2=0')
title('system of 3 equations with 2 unkonwns')

% stepwise solution
figure(2)
Q = zeros(2,2);
F = eye(2,2);
H0 = [2,3]; H1 = [4,3]; H2 = [1,-1]; Bu = 0;
x0 = [0; 0]; P0 = 1000*eye(2,2);
draw_ellipse(x0,P0,0.9,'red',2)
hold on;   
[x1,P1] = kalman(x0,P0,Bu,18,Q,1,F,H0);
draw_ellipse(x1,P1,0.9,'magenta',2)
[x2,P2] = kalman(x1,P1,Bu,24,Q,4,F,H1);
draw_ellipse(x2,P2,0.9,'green',2)
[x3,P3] = kalman(x2,P2,Bu,0,Q,4,F,H2);
draw_ellipse(x3,P3,0.9,'blue',6)

% combined solution
H = [H0; H1; H2]; R = diag([1 4 4]);
z = [18; 24; 0];
[x,P] = kalman(x0,P0,Bu,z,Q,R,F,H);
draw_ellipse(x,P,0.9,'red',2)
axis equal

% plot the iterations for the stepwise solution
figure(3)
hold on;
plot(x0(1),x0(2),'ko');
plot(x1(1),x1(2),'bo');
plot(x2(1),x2(2),'ro');
plot(x3(1),x3(2),'go');
xp = [x0(1),x1(1),x2(1),x3(1)];
yp = [x0(2),x1(2),x2(2),x3(2)];
plot(xp,yp,'k-');
