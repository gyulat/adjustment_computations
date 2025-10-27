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

% number of terms in the mean
nterms = [2,3,4,6,9,25,100,400,2000,10000];
nmax = max(nterms);
% draw enough times, average nterms, estimate the third quartile (75%)
q = 0.75; 
quantile_val = zeros(1,length(nterms));
xi = 0;
k = 1;
x = zeros(1,nsamples);
for i=1:nmax
    x1 = trnd(dof,1,nsamples);
    if i==1
        q1 = quantile(x1, q);
        AE = std(x1);
    end
    x = x + x1;
    if i==nterms(k)
        xav = x/i;
        quantile_val(k) = quantile(xav, q);
        k = k + 1;
    end
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
title(sprintf('Fulfilment of LLN for the %s PDF',name))

grid on
hold on
xl = xlim(); x=0:0.1:xl(2); y=x;
plot(x,y,'color','blue','linestyle','--','linewidth',2)
