% Number of samples to generate
N = 10000;

% Target distribution: standard normal pdf
p = @(x) (1/sqrt(2*pi)) * exp(-0.5*x.^2);

% Proposal distribution: Uniform[-5,5]
a = -5; 
b = 5;
q = @(x) 1/(b-a) * (x >= a & x <= b);  % pdf of uniform
sample_q = @() (a + (b-a)*rand);       % draw from uniform

% Scaling constant M (make sure M*q(x) >= p(x))
% The maximum of p(x) is at x=0 → 1/sqrt(2π) ≈ 0.399
M = 5.0;  % since q(x)=0.1, M=0.399/0.1 ≈ 3.99 < 5 is safe

% Rejection sampling
samples = zeros(N,1);
count = 0;
while count < N
    x_candidate = sample_q();
    u = rand;
    if u <= p(x_candidate)/(M*q(x_candidate))
        count = count + 1;
        samples(count) = x_candidate;
    end
end

% Plot results
figure;
histogram(samples, 50, 'Normalization','pdf');
hold on;
fplot(p, [a b], 'r', 'LineWidth', 2);
legend('Rejection samples','Target N(0,1)');
title('Rejection Sampling Example');
xlabel('x'); ylabel('Density');

