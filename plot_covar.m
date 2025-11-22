function plot_covar(mu,S,col,pconf)
% plot confidence ellipse at level pconf
% mu, S are mean and covariance matrix
% col is the color of the ellipse
% or with 3 arguments only, plot standard ellipse
pkg load statistics  % for Octave
if (nargin==3)
    % standard error ellipse pconf = 0.3935
    pconf = chi2cdf(1,2); 
end
lam = chi2inv(pconf,2);
A = sqrtm(S);
t = [0 : 2*pi/200 : 2*pi];
y = [cos(t); sin(t)];
x = sqrt(lam)*A*y + mu*ones(1,length(t));
plot(x(1,:), x(2,:),col);



