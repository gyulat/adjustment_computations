%% Efficiency of mean and MFV estimations for arbitrary input data vector
clc; clear all; close all

% pkg load statistics % for Octave

% Set random seed for reproducibility
rng(42);

filename = 'madr_spp.txt';  % GPS SPP residuals for MADR station
data = load(filename);
select = data(:,3);
name = 'residual U component';
% select = data(:,1);
% name = 'residual N component';
% select = data(:,2);
% name = 'residual E component';

% Number of samples
nsamples = length(select);

% number of terms in the mean
nterms = [2,3,4,6,9,25,100,400];
nns = length(nterms);
nmax = max(nterms);
q = 0.75; 
x1 = select;  % selected dataset
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

plot(1./sqrt(nterms),quantile_val_mean/q1,'ko','MarkerFaceColor','red')
hold on
plot(1./sqrt(nterms),quantile_val_mfv/q1,'r*')
xlabel('1/√n');
ylabel('q_n/q_1');
title(sprintf('Quantiles from mean and MFV estimation for %s dataset',name))
grid on
legend('mean', 'MFV')

