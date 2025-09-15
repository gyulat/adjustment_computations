function permutation_test(x, y)
    % Permutation test for ratio of standard deviations

    sizex = length(x);
    sizey = length(y);

    pooled = [x; y];
    ratio = std(y) / std(x);

    numSamples = 10000;
    estimates = zeros(numSamples,1);

    % Permutation test
    for k = 1:numSamples
        pooled = pooled(randperm(length(pooled))); % shuffle
        starx = pooled(1:sizex);
        stary = pooled(end-sizey+1:end);
        estimates(k) = std(stary) / std(starx);
    end

    diffCount = sum(estimates <= ratio);
    pvalue = (diffCount / numSamples);

    nbins = ceil(sqrt(numSamples));
    hist(estimates, nbins)
    xlabel('ratio of standard deviations');
    title('Permutation test');
    ylimits = ylim();
    xlimits = xlim();
    text(0.80*xlimits(2), 0.45*ylimits(2), sprintf('sample std ratio: %.3f', ratio));
    text(0.85*xlimits(2), 0.5*ylimits(2), sprintf('p-value: %.3f', pvalue));
    saveas(gcf, ['permutaton-test-result.png']);

    % Print results
    fprintf('permutation test\n');
    fprintf('ratio of std''s data 2./data 1.: %.4f\n', ratio);
    fprintf('number of permutations: %d\n', numSamples);
    fprintf('test result (p-value): %.3f\n', pvalue);
end

