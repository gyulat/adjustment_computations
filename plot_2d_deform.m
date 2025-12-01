%% Plot 2D deformation of a network
clear all; close all

load deform2d.mat

% Used workspace variables:
% X1, Y1, X2, Y2: (npts,1) vectors of coordinates of the network at two epochs
% edg(max npts,2) edges to be plotted, between point indices in rows
%  if not given, plot all edges
% Cov1, Cov2 (2*npts,2*npts) covariance matrices at both epochs
%  1:npts - for X coordinates, npts+1:end - for Y coordinates
% stable(max npts): indices of stable points of the network
% Points{npts}: cell array of point labels

pkg load statistics   % for Octave

function plot_covar(mu,S,col,pconf,filled)
    % plot confidence ellipse at level pconf
    % mu, S are mean and covariance matrix
    % col is the color of the ellipse
    % if filled is true, fill the ellipse
    % or with 3 arguments only, plot standard ellipse
    pkg load statistics  % for Octave
    if (nargin==3)
        % standard error ellipse pconf = 0.3935
        pconf = chi2cdf(1,2); 
    end
    if (nargin==4)
        filled=false;
    end
    lam = chi2inv(pconf,2);
    A = sqrtm(S);
    t = [0 : 2*pi/200 : 2*pi];
    y = [cos(t); sin(t)];
    x = sqrt(lam)*A*y + mu*ones(1,length(t));
    if filled
        fill(x(1,:), x(2,:), col, 'FaceAlpha', 0.7);
    else
        plot(x(1,:), x(2,:), col);
    end
end

pconf = 0.95;  % confidence ellipses

npts = length(X1);  % number of points

% deformations in mm
dX12 = 1000*(X2 - X1);
dY12 = 1000*(Y2 - Y1);

% covariances of deformations
C12 = Cov1 + Cov2;

% plot network edges
figure(1); hold on
if ~exist('edg','var')
    edg = nchoosek(1:npts,2);
end
for i=1:size(edg,1)
    plot([Y1(edg(i,1)),Y1(edg(i,2))],[X1(edg(i,1)),X1(edg(i,2))],'k-')
end

axis equal
xlabel("Y"); 
ylabel("X");

scell = 10;  % scaling of ellipses and deformation vectors


% plot ellipses and deformation vectors
color = 'k';
for i=1:npts
    if any (stable == i)
        filled = true;
    else
        filled = false;
    end
    is = (2*i):-1:(2*i-1);
    plot_covar([Y1(i);X1(i)],scell^2*C12(is,is),color,pconf,filled)
    text(Y1(i)+3,X1(i)+3,Points{i})
    quiver(Y1(i),X1(i),scell*dY12(i),scell*dX12(i),0,'red')
end
scatter(Y1,X1,10,'k','filled','o')


xlimits = xlim(); dx = diff(xlimits);
ylimits = ylim(); dy = diff(ylimits);

% plot scale
plot([xlimits(2)-0.9*dx,xlimits(2)-0.9*dx+scell],[ylimits(2)-0.1*dy,ylimits(2)-0.1*dy],'k-')
text(xlimits(2)-0.9*dx+0.4*scell,ylimits(2)-0.1*dy+ 0.2*scell, '1 mm')

title(sprintf("%.1f %% confidence ellipses",100*pconf))
