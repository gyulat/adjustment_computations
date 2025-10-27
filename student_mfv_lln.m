%% Law of large numbers for the Student supermodel
clc; clear all; close all

% pkg load statistics % for Octave

% Set random seed for reproducibility
rng(42);

% Number of samples
nsamples = 10000;

% degree of freedom of the Student distribution
dof = 1;
name = 'Cauchy';
dof = 4;
name = 'Jeffreys a=5';

% number of terms in the MFV
nterms = [2,3,4,6,9,25,100,400];
nmax = max(nterms);
% draw enough times, MFV of nterms, estimate the third quartile (75%)
q = 0.75; 
nns = length(nterms); % number of cases
x = trnd(dof,nsamples,nmax);
q1 = quantile(x(:,1),q);
AE = std(x(:,1));
quantile_val = zeros(1,nns);
for i=1:nns
    mfv = MFV_matrix(x(:,1:nterms(i)));
    quantile_val(i) = quantile(mfv, q);
end

% Plot scaled quantiles in terms of 1/sqrt(n)
if dof == 1
    plot(1./sqrt(nterms),quantile_val/q1,'ko','MarkerFaceColor','red')
    ylabel('q_n/q_1');
else
    plot(1./sqrt(nterms),quantile_val/AE/0.6745,'ko','MarkerFaceColor','red')
    ylabel('q_n/(A_E*0.6745)');
end
xlabel('1/√n');
title(sprintf('Fulfilment of LLN for the %s PDF with MFV',name))

grid on
hold on
xl = xlim(); x=0:0.1:xl(2); y=x;
plot(x,y,'color','blue','linestyle','--','linewidth',2)
