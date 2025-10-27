%% Efficiency of mean and MFV estimations for Student data distribution
clc; clear all; close all

% pkg load statistics % for Octae

% Set random seed for reproducibility
rng(42);

% Number of samples
nsamples = 10000;

% degree of freedom of the Student distribution
dof = 1;
name = 'Cauchy';
% dof = 4;
% name = 'Jeffreys a=5';

% number of terms in the mean
nterms = [2,3,4,6,9,25,100,400];
nns = length(nterms);
nmax = max(nterms);
q = 0.75; 
x1 = trnd(dof,1,nsamples);  % draw one sample
% AE = std(x1);
q1 = quantile(x1, q);

% resample enough times, mean/MFV of nterms, estimate the third quantiles (75%)
xav = zeros(nsamples,nns);
xmfv = zeros(nsamples,nns);
quantile_val_mean = zeros(1,nns);
quantile_val_mfv = zeros(1,nns);
for i=1:nsamples
    xr = randsample(x1,nsamples,true);  % resample
    for k=1:nns
        xav(i,k) = mean(xr(1:nterms(k)));
        xmfv(i,k) = MFV1(xr(1:nterms(k)));
    end
end
for k=1:nns
    quantile_val_mean(k) = quantile(xav(:,k), q);
    quantile_val_mfv(k) = quantile(xmfv(:,k), q);
end

% Plot third quartiles in terms of 1/sqrt(n)
figure(1);
plot(1./sqrt(nterms),quantile_val_mean/q1,'ko','MarkerFaceColor','red')
hold on
plot(1./sqrt(nterms),quantile_val_mfv/q1,'r*')
xlabel('1/√n');
ylabel('q_n/q_1');
title(sprintf('Quantiles from mean and MFV estimation for %s PDF',name))
grid on
legend('mean', 'MFV')

% Plot quantile ratios of MFV/mean
figure(2)
plot(1./sqrt(nterms),quantile_val_mfv./quantile_val_mean,'b*')
xlabel('1/√n');
ylabel('q_{MFV}/q_{mean}');
title(sprintf('Quantile ratios of MFV and mean estimation for %s PDF',name))
grid on

% Plot relative efficiency of MFV/mean
figure(3)
plot(1./sqrt(nterms),(quantile_val_mfv./quantile_val_mean).^2,'b*')
xlabel('1/√n');
ylabel('(q_{MFV}/q_{mean})^2');
title(sprintf('Relative efficiency of mean estimation wrt. to MFV for %s PDF',name))
grid on
