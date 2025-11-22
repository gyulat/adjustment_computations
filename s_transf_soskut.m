%% S-transformation of the solution
%  plot network and error ellipses
clear all; close all
load soskut.mat

pconf = 0.3935; % standard error ellipses
pconf = 0.95;  % confidence ellipses

npts = length(XYZ);  % number of points

% fix 3 coordinates: Point 6: x, y; Point 4: x 
% indexes of fixed coordinates
idx = [Indexes(4,1),Indexes(6,1),Indexes(6,2)];
% datum selector matrix E
d = zeros(1,12); d(idx) = 1;
E = diag(d);

% set up G-matrix
G = zeros(12,3);
% separate coordinates
X = XYZ(:,1); Y = XYZ(:,2);
% indexes of points
ix = Indexes(:,1); iy = Indexes(:,2); 
ip = fix(Indexes(:,2)/2);
G(ix,1) = 1; G(iy,2) = 1;
G(ix,3) = -X(ip); G(iy,3) = Y(ip);

% S-matrix
S = eye(12) - G*inv(G'*E*G)*G'*E;

% solution
xS = S*x(1:12);

% adjusted coordinates
XS = XYZ_0(:,1)+0.001*xS(ix);
YS = XYZ_0(:,2)+0.001*xS(iy);
fprintf('Point     X          Y \n')
for i=1:6
    fprintf('  %d  %10.4f %10.4f\n',i,XS(i),YS(i));
end

%    Point     X          Y 
%      1   1372.4832   564.3009
%      2   1370.3007   404.0361
%      3   1326.4009   188.2275
%      4   1000.0000   607.0909
%      5    974.9404   383.4119
%      6   1000.0000     0.0000


% covariance
CS = S*C_xx(1:12,1:12)*S';

% plot network
figure(1); hold on
edg = nchoosek(1:npts,2); 
for i=1:length(edg)
    plot([Y(edg(i,1)),Y(edg(i,2))],[X(edg(i,1)),X(edg(i,2))],'k-')
end
axis equal
xlabel("Y"); 
ylabel("X");
xlimits = xlim(); dx = diff(xlimits);
ylimits = ylim(); dy = diff(ylimits);
scell = 100;  % scaling up ellipse

% plot scale
plot([xlimits(2)-0.9*dx,xlimits(2)-0.9*dx+scell],[ylimits(2)-0.1*dy,ylimits(2)-0.1*dy],'k-')
text(xlimits(2)-0.9*dx+0.4*scell,ylimits(2)-0.1*dy+ 0.2*scell, '1 mm')

colors = ['r','g','b','m','k','c'];

for i=1:npts
    i1 = Indexes(i,1);i2 = Indexes(i,2); is = i2:-1:i1;
    plot_covar([Y(i);X(i)],scell^2*CS(is,is),colors(i),pconf)
    text(Y(i)+10,X(i)-10,Points{i})
end
title(sprintf("%.1f %% confidence ellipses",100*pconf))
