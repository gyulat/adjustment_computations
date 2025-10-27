%% Law of large numbers for the Tukey supermodel PDF with bootstrap sampling
clc; clear all; close all

% pkg load statistics  % for Octave

% Set random seed for reproducibility
rng(42);

% Number of samples
nsamples = 10000;

% parameters of the Tukey distribution
p = 0.25;
% p = 0.0;  % for a Gaussian
sc = 150;

% number of terms in the mean
nterms = [2,3,4,6,9,25,100,400,2000,10000];
nns = length(nterms);
nmax = max(nterms);
x1 = draw_tukey(p,sc,nsamples); % draw one sample
AE = std(x1);
% resample enough times, average nterms, estimate the third quantile (75%)
q = 0.75; 
xav = zeros(nsamples,nns);
quantile_val = zeros(1,nns);
for i=1:nsamples
    xr = randsample(x1,nsamples,true);  % resample
    for k=1:nns
        xav(i,k) = mean(xr(1:nterms(k)));       
    end
end
for k=1:nns
    quantile_val(k) = quantile(xav(:,k), q);
end

% Plot scaled quantiles in terms of 1/sqrt(n)

plot(1./sqrt(nterms),quantile_val/0.6745/AE,'ko','MarkerFaceColor','red')
xlabel('1/√n');
ylabel('q_n/(A_E*0.6745)');
if abs(p-0.25) < 1e-6
    title('Fulfilment of LLN for the Tukey supermodel (p=0.25, σ_c=150 with bootstrap)')
else
    title('Fulfilment of LLN for the standard Gaussian with bootstrap')
end
grid on
hold on
xl = xlim(); x=0:0.1:xl(2); y=x;
plot(x,y,'color','blue','linestyle','--','linewidth',2)
