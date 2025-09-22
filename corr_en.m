% Load data, skipping the first row
data = dlmread('values.txt', '', 1, 0);

% Number of data points (rows) and variables (columns)
[nd, nv] = size(data);

% Covariance matrix
covmx = cov(data);

% Variable labels
v = {'x', 'y'};

% Display results
fprintf('Covariance matrix from Monte Carlo simulation\n');
fprintf(' number of data: %d\n', nd);
fprintf('covariance matrix:\n');

for i = 1:nv
    for j = 1:nv
        fprintf('%2s%2s: %8.4e\n', v{i}, v{j}, covmx(i,j));
    end
end




