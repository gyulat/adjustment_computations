pkg load statistics % Octave - for Matlab use Statistics Toolbox
%% Example 1: test for mean value
clc; clear all; close all
fprintf('\nTest for mean value\n\n');
% measured angles are 180-00-00 + da, where
da = [4.8    10.0   -6.4    6.0    3.9];  % arcsec
damean = mean(da); % 3.6600
% u-statistics for equality test with 180-00-00 and std 4 arcsec
astd = 4;
u = -damean/(astd/sqrt(length(da)));  % -2.0460
% confidence level (two-sided test)
pconf = 0.90;
% test statistic from inverse normal distribution
ut = norminv(pconf+0.5*(1-pconf)); % 1.6449
if (abs(u) < ut)
    disp('Fail to reject the null hypothesis: mean is equal to 180-00-00');
else
    disp('Reject the null hypothesis: mean is not equal to 180-00-00');
end
fprintf('Mean value: 180-00-%04.1f\n', damean);
% Display the calculated values for reference
fprintf('Test Statistic (u): %.4f\n', u);
fprintf('Critical Value (ut): %.4f\n', ut);
% Reject the null hypothesis: mean is not equal to 180-00-00
% test for mean value with empirical standard deviation
estd = std(da); % 6.0875
t = -damean/(estd/sqrt(length(da))); % -1.3444
% critical value from t-distribution for two-sided test
df = length(da) - 1; % degrees of freedom
tt = tinv(pconf + 0.5*(1 - pconf), df); % critical t-value
if (abs(t) < tt)
    disp('Fail to reject the null hypothesis: mean is equal to 180-00-00 with empirical std');
else
    disp('Reject the null hypothesis: mean is not equal to 180-00-00 with empirical std');
end
% Display the calculated values for reference
fprintf('Test Statistic (t): %.4f\n', t);
fprintf('Critical Value (tt): %.4f\n', tt);
pause

%% Example 2: deformation measurement
% deformation detection measurements for dam movement
clear all; close all
fprintf('\nDeformation measurement statistical test\n\n');
% epoch 1, coordinate in the direction perpendicular to the dam
a1 = 100.4136;
% number of measurements
k1 = 7;
% empirical std
s1 = 0.0042;
% epoch 2, coordinate in the direction perpendicular to the dam
a2 = 100.4184;
% number of measurements
k2 = 5;
% empirical std
s2 = 0.0036;
% Calculate the difference in measurements between epochs
deltaA = a2 - a1;
% test for equality of two mean values
% Calculate the pooled standard deviation
sp = sqrt(((k1 - 1) * s1^2 + (k2 - 1) * s2^2) / (k1 + k2 - 2));
% Calculate the t-statistic for the difference in means
tDiff = deltaA / (sp * sqrt(1/k1 + 1/k2));
pconf = 0.90;
fprintf('Test with confidence level %.2f\n\n',pconf);
% critical value from t-distribution for two-sided test
ttDiff = tinv(pconf + 0.5*(1 - pconf), k1 + k2 - 2); % critical t-value for difference
if (abs(tDiff) < ttDiff)
    disp('Fail to reject the null hypothesis: means are equal between epochs');
else
    disp('Reject the null hypothesis: means are not equal between epochs');
end
% Display the calculated values for reference
fprintf('Pooled Standard Deviation: %.4f\n', sp);
fprintf('Difference in Measurements: %.4f\n', deltaA);
fprintf('t-Statistic for Difference: %.4f\n', tDiff);
fprintf('Critical t-Value for Difference: %.4f\n', ttDiff);
% Calculate the p-value for the t-test
pValue = 2 * (1 - tcdf(abs(tDiff), k1 + k2 - 2));
fprintf('P-Value for Difference: %.4f\n', pValue);
% Display the conclusion based on the p-value
if pValue < 0.05
    disp('Reject the null hypothesis: means are significantly different.');
else
    disp('Fail to reject the null hypothesis: means are not significantly different.');
end

% repeat the test with higher confidence level
pconf = 0.95; 
fprintf('\nRepeat the test with higher confidence level %.2f\n\n',pconf);
% critical value from t-distribution for two-sided test
ttDiff = tinv(pconf + 0.5*(1 - pconf), k1 + k2 - 2); % critical t-value for difference
if (abs(tDiff) < ttDiff)
    disp('Fail to reject the null hypothesis: means are equal between epochs');
else
    disp('Reject the null hypothesis: means are not equal between epochs');
end
fprintf('Critical t-Value for Difference: %.4f\n', ttDiff);
pause



%% Example 3: test for empirical standard deviation
clear all; close all
fprintf('\nStatistical test for empirical standard deviation\n\n');
% Calculate the empirical standard deviation for the next test
da = [4.8    10.0   -6.4    6.0    3.9];  % arcsec
sEmpirical = std(da); % Calculate empirical standard deviation
fprintf('Empirical Standard Deviation: %.4f\n', sEmpirical);
sTheoretical = 4.0;
% Calculate the chi-squared statistic for the empirical standard deviation
chiSquared = (length(da) - 1) * (sEmpirical^2 / sTheoretical^2);
% Calculate the critical values for the chi-squared distribution
pconf = 0.95;
fprintf('Test with confidence level %.2f\n\n',pconf);
chiCriticalLow = chi2inv((1-pconf)/2, length(da) - 1);
chiCriticalHigh = chi2inv(pconf/2, length(da) - 1);
% Determine if the empirical standard deviation is significantly different from the theoretical value
if chiSquared < chiCriticalLow || chiSquared > chiCriticalHigh
    disp('Reject the null hypothesis: empirical std is significantly different from theoretical std.');
else
    disp('Fail to reject the null hypothesis: empirical std is not significantly different from theoretical std.');
end
fprintf('Chi-Squared Statistic: %.4f\n', chiSquared);
fprintf('Critical Values: [%.4f, %.4f]\n', chiCriticalLow, chiCriticalHigh);

fprintf('\nEnd of tests\n');

