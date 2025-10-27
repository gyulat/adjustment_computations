%% Law of large numbers for the Tukey supermodel PDF
clc; clear all; close all

%pkg load statistics  % for Octave

% Set random seed for reproducibility
rng(42);

% Number of samples
nsamples = 10000;

% parameters of the Tukey supermodel PDF (from Steiner)
p = 0.25;
% p = 0; % for Gaussian PDF
sc = 150;

% number of terms in the mean
nterms = [2,3,4,6,9,25,100,400,2000,10000];
nmax = max(nterms);
% draw enough times, average nterms, estimate the third quantile (75%)
q = 0.75; 
quantile_val = zeros(1,length(nterms));
xi = 0;
k = 1;
x = zeros(1,nsamples);
for i=1:nmax
    x1 = draw_tukey(p,sc,nsamples);
    if i==1
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

plot(1./sqrt(nterms),quantile_val/0.6745/AE,'ko','MarkerFaceColor','red');
xlabel('1/√n');
ylabel('q_n/(A_E*0.6745)');
if abs(p-0.25) < 1e-6
    title('Fulfilment of LLN for the Tukey supermodel (p=0.25, σ_c=150)')
else
    title('Fulfilment of LLN for the standard Gaussian')
end
grid on
hold on
xl = xlim(); x=0:0.1:xl(2); y=x;
plot(x,y,'color','blue','linestyle','--','linewidth',2)
