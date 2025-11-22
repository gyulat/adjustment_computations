%% plot network and error ellipses
clear all; close all
load pillar4.mat

% separate coordinates
x = XYZ(:,1); y = XYZ(:,2);
npts = length(x);  % number of points

% --- set confidence level
pconf = 0.3935; % standard error ellipses
pconf = 0.99;  % confidence ellipses

% plot network - edges in all combinations
figure(1); hold on
edg = nchoosek(1:npts,2); 
for i=1:length(edg)
    plot([y(edg(i,1)),y(edg(i,2))],[x(edg(i,1)),x(edg(i,2))],'k-')
end
axis equal
xlabel("Y"); 
ylabel("X");
xlimits = xlim(); dx = diff(xlimits);
ylimits = ylim(); dy = diff(ylimits);

% --- set this scale according to your needs
scell = 10;  % scaling up the ellipse

% --- plot scale - customize if necessary
plot([xlimits(2)-0.9*dx,xlimits(2)-0.9*dx+scell],[ylimits(2)-0.1*dy,ylimits(2)-0.1*dy],'k-')
text(xlimits(2)-0.9*dx+0.4*scell,ylimits(2)-0.1*dy+ 0.2*scell, '1 mm')

colors = ['r','g','b','m','k','c'];

for i=1:npts
    i1 = Indexes(i,1);i2 = Indexes(i,2); is = i2:-1:i1;
    plot_covar([y(i);x(i)],scell^2*C_xx(is,is),colors(i),pconf)
    text(y(i)+1.0,x(i)-1.0,Points{i})
end
title(sprintf("%.1f %% confidence ellipses",100*pconf))


