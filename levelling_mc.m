% levelling_mc.m
% Calculate standard deviation of heights in a levelling network by Monte Carlo simulation
%
% Network consists of points 3, 6, 9. 
% Height of point 3 is 0.00000 m, fixed

% height differences from adjusted point heights are the measurements
obs = [7.44515, -8.31910, 0.87420];  % m
% standard errors for these measurements
% distances levelled
d = [580.5, 445.3, 511.2]; % m
% standard error of levelling per 1 km double-run
m0 = 0.4;  % mm/km
ost = m0*sqrt(d/1000);  % 0.3048, 0.2669, 0.2860

% known misclosure of the adjusted levelling network
misclosure = 0.4;  % mm

% Monte Carlo experiment
rng(0);               % for reproducibility
Nsim = 200000;        % number of simulated set of measurements
% the simulated height differences have a Gaussian PDF with known mean and std and sum to misclosure for each trial
n = length(obs);               % number of variables
mu = obs';  % mean
sigma = ost/1000;  % standard deviations in m
variances = sigma.^2;
Vsum = sum(variances);
c = misclosure/1000;  % fixed value of misclosure in m

% === Simulation: draw X and apply adjustment ===
% Draw X: shape (n x Nsim)
X = repmat(mu,1,Nsim) + sigma' .* randn(n, Nsim);

% Compute adjustment factor for each draw: scalar per sample
adjustment = (c - sum(X,1)) / Vsum;

% Form Y: Y = X + variances * adjustment  (variances is n x 1, adjustment is 1 x Nsim)
Y = X + variances' * adjustment;   % results in n x Nsim

% === Empirical estimates ===
emp_mean = mean(Y,2);              % n x 1
emp_cov = cov(Y');                 % n x n   (observations are rows -> Y' is Nsim x n)
% emp_cov

% === Checks and comparisons ===
max_sum_error = max(abs(sum(Y,1) - c));   % should be ~0 (numerical precision)
mean_diff_inf = max(abs(emp_mean - mu));   % max absolute diff in mean
std_diff_norm = norm(sqrt(diag(emp_cov)) - sigma'');  % norm of stdev diff.

% Print results
fprintf('Simulation count: %d\n', Nsim);
fprintf('Max absolute deviation of sample sums from c: %.3e\n', max_sum_error);
fprintf('Max abs difference between empirical mean and prescribed mean: %.3e\n', mean_diff_inf);
fprintf('norm of (empirical std - given std): %.3e\n', std_diff_norm);
fprintf('\n3 components: (empirical_mean, prescribed mean)\n');
disp([emp_mean(1:3), mu(1:3)]);

% calculate adjusted heights from simulated observations
mc = sum(Y,1)-c;
% 1000*mc(1:10) % first 10 misclosures in mm
% weights for measurements are proportional to distances
w = d/sum(d); %  0.3777   0.2897   0.3326
% corrections to observations
cobs = -c*w; 
% corrected simulated observations
Yc = cobs + Y';  % Nsim x n
% check sum of corrected measurements whether it is 0
ss = sum(Yc,2);
% ss(1:3)  %  -4.4409e-16, ... etc - OK
% calculate heights of points 6 and 9
m6 = Yc(:,1); 
m9 = -Yc(:,3);

% calculate mean heights and their std from MC experiment
mm6 = mean(m6);
std6 = std(m6);

mm9 = mean(m9);
std9 = std(m9);
fprintf('Simulation results\n');
fprintf('Point 6: height: %.5f m,   std: %.4f mm\n', mm6, 1000*std6);
fprintf('Point 9: height: %.5f m,   std: %.4f mm\n', mm9, 1000*std9);

% Histograms
figure(1)
hist(1000*(m6-mm6), 50)
title('Benchmark 6')
xlabel('height residuals (mm)')
figure(2)
hist(1000*(m9-mm9), 50)
title('Benchmark 9')
xlabel('height residuals (mm)')


%Simulation count: 200000
%Max absolute deviation of sample sums from c: 1.376e-15
%Max abs difference between empirical mean and prescribed mean: 5.616e-05
%norm of (empirical std - given std): 1.673e-04

%3 components: (empirical_mean, prescribed mean)
%   7.4452   7.4451
%  -8.3191  -8.3191
%   0.8743   0.8742
%Simulation results
%Point 6: height: 7.44506 m,   std: 0.2401 mm
%Point 9: height: -0.87412 m,   std: 0.2330 mm


