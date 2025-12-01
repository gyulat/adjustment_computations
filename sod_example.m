%% SOD example of three measured distances 

% design matrix
A = -[sind(30), cosd(30); sind(140), cosd(140); sind(250), cosd(250)]
sq = A.^2;
pr = A(:,1).*A(:,2);

% set up the equations for unknown weights
S = [sq(:,1)';pr'; sq(:,2)']
b = [1 0 1]';

% solution
p = S\b
%    p =
%       0.5662
%       0.8675
%       0.5662
% check
norm(S*p-b) % 2.2204e-16

% accuracies of single measurements
d = [13.58, 9.15, 6.94]; % km, distances
ss = 0.5 + 0.1*d  % 1.8580   1.4150   1.1940

% required number of measurements
n = p'.*ss.^2  % 1.9547   1.7370   0.8072

% realizable number of measurements
nreal = [2, 2, 1];

% realized accuracies
s = ss./sqrt(nreal)  % 1.3138   1.0006   1.1940

% realized weights
preal = 1./s.^2  % 0.5793   0.9989   0.7014

P = diag(preal);
N = A'*P*A
%    N =
%       1.176942  -0.015553
%      -0.015553   1.102735

% realized cofactor matrix
Q = inv(N)
%    Q =
%       0.849818   0.011986
%       0.011986   0.907005

% plot confidence/standard error ellipse
pkg load statistics   % for Octave

function plot_covar(mu,S,col,pconf,filled)
    % plot confidence ellipse at level pconf
    % mu, S are mean and covariance matrix
    % col is the color of the ellipse
    % if filled is true, fill the ellipse
    % or with 3 arguments only, plot standard ellipse
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

function [a,b,alpha] = std_ellipse(C)
    [s,lam] = eig(C);
    fac = 1;
    alpha = 90/pi*atan(2*C(1,2)/(C(1,1)-C(2,2)));
    if alpha<0
        alpha = alpha + 180;
    end
    a = fac*sqrt(lam(2,2));
    b = fac*sqrt(lam(1,1));
end

[a,b,alpha] = std_ellipse(Q)
%    a = 0.9536
%    b = 0.9205
%    alpha = 168.63

plot_covar([0;0],Q,'r',0.3935)
hold on
plot_covar([0;0],eye(2),'b',0.3935)
axis equal

