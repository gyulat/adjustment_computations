% Simulate standard Gaussian random variables with their sum equals a fixed value
clc; clear all; close all

% Parameters
n = 5;           % number of variables
c = 50;          % desired sum
mu = linspace(0,20,5);         % mean of variables
sigma = 1;      % std of variables

% Ensure mu is column vector
if isscalar(mu)
    mu = mu * ones(n,1);
else
    mu = mu(:);
end

Nsim = 500000;
results = zeros(Nsim,n);
X = zeros(Nsim,n);
Y = zeros(Nsim,n);
ones_vec = ones(n,1);

Vsum = n*sigma;
variances = sigma.^2; 

X = repmat(mu,1,Nsim) + sigma .* randn(n, Nsim);

% Compute adjustment factor for each draw: scalar per sample
% adjustment = (c - sum(X,1)) / Vsum  -> 1 x Nsim
adjustment = (c - sum(X,1)) / Vsum;

% Form Y: Y = X + variances * adjustment  (variances is n x 1, adjustment is 1 x Nsim)
Y = X + variances * adjustment;   % results in n x Nsim

% === Check ===
max_sum_error = max(abs(sum(Y,1) - c));   % should be ~0 (numerical precision)

% Print results
fprintf('Simulation count: %d\n', Nsim);
fprintf('Max absolute deviation of sample sums from c: %.3e\n', max_sum_error);

%    Max absolute deviation of sample sums from c: 1.421e-14

% plot histograms
hold on;
colors = ['g','r','b','m','y'];
for i=1:n
    hist(Y(i,:), 30, 1, 'facecolor', colors(i));
end
hold off;
title(sprintf('%d Gaussians with sum %.0f (%d runs)',n,c,Nsim))
