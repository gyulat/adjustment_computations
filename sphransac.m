function sphransac
    % RANSAC sphere fitting
    % clc; clear all
    % test data
    X = load('sphere.dat');
    size(X)
    nd = length(X);
    tol = 0.003;  % threshold for fit
    np = 4; % number of parameters
    k = 25;  % number of iterations
    nmax = 0; % empty consensus set
    for i=1:k
      % select np points at random
      is = randperm(nd,np);
      % determine parameters of the sphere
      p = sph4fit(X(is,:));
      if length(p)==0  % no solution
        continue
      end
      % data distances to sphere
      t = sphdist(p, X);
      ins = t<tol;
      Xk = X(ins);  % conform data
      nin = length(Xk);  % cardinality of the consensus set  
      if nin > nmax  % so far the best
        Xin = Xk; nmax = nin; inliers = ins;
        bp = p;  % best sphere
      end
    end
    fprintf("number of inliers : %d\n",nmax);
    fprintf("number of outliers: %d\n",nd-nmax);
    p = sph4fit(X(inliers,:));  % LSQ best fitting sphere
    fprintf("parameters of the best fitting sphere :\n");
    fprintf('centre x,y,z      : %f, %f, %f\n', p(1),p(2),p(3));
    fprintf('radius            : %f\n', p(4));
    % plot results
    figure(1);
    size(X(inliers,:))
    size(X(~inliers,:))
    plot3(X(inliers,1),X(inliers,2),X(inliers,3),'g.')
    hold on
    plot3(X(~inliers,1),X(~inliers,2),X(~inliers,3),'r.')
    axis equal
    % Plot histogram of residuals of conform points
    % calculate signed distances (residuals)
    v = sphdist(p,X(inliers,:),0);
    figure(2);
    hist(v, 50);
    % check normality with a quantile-quantile plot
    %pkg load statistics  % for Octave
    figure(3);
    qqplot(v);
end


function x = svdsolve(A, b)
    % Solve linear system using SVD
    [U,S,V] = svd(A);
    s = diag(S);
    c = U' * b;
    w = c(1:length(s)) ./ s;
    x = V * w;
end


function p = sph4fit(X)
    % sph4fit(X) fits a sphere to four or more 3D points.
    % Output: [x0, y0, z0, r]

    % Subtract mean
    Xm = mean(X,1);
    Xs = X - Xm;
    n = size(Xs,1);

    % Construct matrix M and vector h
    M = [-2*Xs, ones(n,1)];
    h = -sum(Xs.^2, 2);

    % Check conditioning
    if cond(M) > 1000
        p = [];
        return;
    end

    % Solve M*p = h
    sol = svdsolve(M, h);

    % Radius and centre
    r = sqrt(sum(sol(1:3).^2) - sol(4));
    c = sol(1:3)' + Xm;

    % Output vector
    p = [c, r];
end


function dist = sphdist(sph, X, pos)
    % sphdist(sph,X,pos): distance of points X from sphere sph
    if nargin < 3
        pos = 1;
    end
    r = sph(4);
    c = sph(1:3);

    % Distance from centre squared
    dc = sum((X - c).^2, 2);

    if pos == 1
        dist = abs(sqrt(dc) - r);
    else
        dist = sqrt(dc) - r;
    end
end
