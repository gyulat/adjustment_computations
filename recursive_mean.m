% Kalman filter: recursive estimation of mean
clc; clear all; close all
npt = 200;
x = 5 + randn(1,npt);   % measurements
n = 1:npt;              % time steps
K = 1 ./ (n+1);         % Kalman gain
mu = zeros(1,npt);      % estimated mean
mu(1) = x(1);

for i = 2:npt
    mu(i) = mu(i-1) + K(i-1) * (x(i) - mu(i-1));
end

plot(n, x, 'o', 'DisplayName', 'measurements'); hold on;
plot(n, mu, 'LineWidth', 2, 'DisplayName', 'estimation of mean');
legend show;
title('Recursive estimation of mean');
hold off;

