function cylransac
    % RANSAC cylinder fitting
    % clc; clear all
    % test data
    X = load('cyl.dat');
    nd = length(X);
    tol = 10;  % threshold for fit
    np = 5; % number of parameters
    k = 25;  % number of iterations
    nmax = 0; % empty consensus set
    for i=1:k
      % select np points at random
      is = randperm(nd,np);
      % determine parameters of the cylinder
      p = cyl5fit(X(is,:));
      % number of cylinders
      if p==0
        continue
      else
        nsol = size(p,2);
      end
      % fprintf("iteration %d, %d solutions\n",i,nsol);
      % pause()
      for j=1:nsol
          % data distances to cylinder
          t = cyldist(p(:,j), X);
          ins = t<tol;
          Xk = X(ins);  % conform data
          nin = length(Xk);  % cardinality of the consensus set  
          if nin > nmax  % so far the best
            Xin = Xk; nmax = nin; inliers = ins;
            bp = p;  % best cylinder
          end
      end
    end
    fprintf('number of inliers    : %d\n',nmax);
    fprintf('number of outliers   : %d\n',nd-nmax);
    [c, r, a] = cylinder(X(inliers,:));
    fprintf('point on the axis    : [ %f %f %f ] \n', c);
    fprintf('radius               : %f\n', r);
    fprintf('axis (theta, lambda) : [ %f %f ] \n', a(1),a(2));
    % plot results
    figure(1);
    plot3(X(inliers,1),X(inliers,2),X(inliers,3),'g.')
    hold on
    plot3(X(~inliers,1),X(~inliers,2),X(~inliers,3),'r.')
    axis equal
end


function cyl = cyl5fit(X)
    % fits cylinders to 5 points stored as rows of matrix X
    % Reference: Paul (2006): A Cylinder of Revolution on Five Points
    % Results: parameters stored in a vector of 7 components  
    %            cyl: [r,ax,ay,az,fx,fy,fz], where:
    %              r: radius of the cylinder
    %              ax,ay,az: unit vector of the cylinder's axis
    %              fx,fy,fz: coordinates of a point on the axis
    %    Remark:  When there is no good solution, the function returns False

    % Shift to origin
    o=X(1,:)-X(1,:);
    pt=X(2,:)-X(1,:);
    qt=X(3,:)-X(1,:);
    rt=X(4,:)-X(1,:);
    st=X(5,:)-X(1,:);

    % Rotation such that OP=x and OQ lies in xy-plane
    x = pt / norm(pt);
    y = qt - (dot(qt, pt) / norm(pt)) * x;
    y = y / norm(y);
    z = cross(x, y);
    R = [x; y; z];

    p = R * pt';
    q = R * qt';
    r = R * rt';
    s = R * st';

    p1 = p(1); q1 = q(1); q2 = q(2);
    r1 = r(1); r2 = r(2); r3 = r(3);
    s1 = s(1); s2 = s(2); s3 = s(3);

    % --- Equation 11 ---
    b1 = 2*q2*r2*(q1 - r1);
    b2 = q2*r3*(p1 - 2*q1);
    b3 = q2*r3*(p1 - 2*r1);
    b4 = q2*(r1*(r1 - p1) + r2*(r2 - q2)) + q1*r2*(p1 - q1);
    b5 = r3*(q2*(q2 - 2*r2) + q1*(q1 - p1));
    b6 = q2^2 * r3;
    b7 = q1*r3*(q1 - p1);
    b8 = q2*(r2*(r2 - q2) + r3^2);
    b9 = q2*(r1*(r1 - p1) + r3^2) + q1*r2*(p1 - q1);

    c1 = q2*(r2*s3*(r2 - q2) + r3*s2*(q2 - s2) + r3*s3*(r3 - s3));
    c2 = q2*r3*s3*(r3 - s3) + q1^2*(r3*s2 - r2*s3) + ...
         p1*(q2*(r3*s1 - r1*s3) + q1*(r2*s3 - r3*s2)) + ...
         q2*(r1^2*s3 - r3*s1^2);
    c3 = q2*(p1*(s1*r3 - s3*r1) - r3*(s1^2 + s2^2) + s3*(r1^2 + r2^2)) + ...
         (r2*s3 - r3*s2)*(q1*(p1 - q1) - q2^2);
    c4 = 2*q2*r3*s3*(s2 - r2);
    c5 = 2*q2*r3*s3*(s1 - r1);
    c6 = 2*q2*(r3*s2*(s1 - q1) + r2*s3*(q1 - r1));

    % --- Equation 10 ---
    d1 = c1*b2 - b6*c6;
    d2 = c6*b7 - b2*c2;
    d3 = c1*b7 - b6*c2;
    d4 = c1*b1 - b8*c6 - b6*c5;
    d5 = c5*b7 + c6*b9 - b1*c2 - b2*c4;
    d6 = c5*b9 + c6*b5 - b3*c2 - b1*c4 - b2*c3;
    d7 = c1*b5 - b8*c4 - b6*c3;
    d8 = c1*b9 - b8*c2 - b6*c4;
    d9 = c1*b3 - b8*c5;
    d10 = c5*b5 + c6*b4 - b3*c4 - b1*c3;
    d11 = c1*b4 - b8*c3;
    d12 = c5*b4 - b3*c3;

    % --- Equation 9 ---
    D1 = d1*d2 - d3^2;
    D2 = d4*d2 + d1*d5 - 2*d8*d3;
    D3 = d9*d2 + d4*d5 + d1*d6 - 2*d7*d3 - d8^2;
    D4 = d9*d5 + d4*d6 + d1*d10 - 2*d11*d3 - 2*d7*d8;
    D5 = d9*d6 + d4*d10 + d1*d12 - 2*d11*d8 - d7^2;
    D6 = d9*d10 + d4*d12 - 2*d11*d7;
    D7 = d9*d12 - d11^2;

    % --- Equation 8 ---
    pol6 = [D1 D2 D3 D4 D5 D6 D7];
    r6 = roots(pol6);

    % Real roots only
    a2_all = real(r6(abs(imag(r6)) < 1e-12));
    n = numel(a2_all);
    if n == 0
        cyl = 0;
        return
    end

    % Ensure column vectors for a1,a2 calculations
    a2 = a2_all(:);           % n x 1 column
    % --- Compute a1 (as column) ---
    e1 = c1*(b2*a2.^2 + b1*a2 + b3) - (b6*a2 + b8).*(c6*a2 + c5);
    e2 = c1*(b7*a2.^3 + b9*a2.^2 + b5*a2 + b4) - ...
         (b6*a2 + b8).*(c2*a2.^2 + c4*a2 + c3);
    a1 = (-e2 ./ e1);        % n x 1 column

    % Build 'a' as 3 x n (rows: a1, a2, 1)
    a = [a1.'; a2.'; ones(1, n)];   % 3 x n
    norma = sqrt(sum(a.^2, 1));     % 1 x n

    % --- Equation (5): compute f (keep f1,f2,f3 as column vectors then stack) ---
    asq = a1.^2 + a2.^2 + 1;        % n x 1
    f1 = p1*(1 + a2.^2) ./ (2*asq); % n x 1
    f2 = (q1^2 + q2^2 + (q1.*a2 - q2.*a1).^2) ./ (2*asq); % n x 1
    f2 = (f2 - q1.*f1) / q2;
    f3 = -(a1.*f1 + a2.*f2);       % n x 1

    % Stack into 3 x n matrix (transpose each column vector)
    f = [f1.'; f2.'; f3.'];        % 3 x n   <-- <== your suggested fix

    % Cylinder radius (1 x n)
    radius = sqrt(sum(f.^2, 1));

    % Back-transform (rotation and shift)
    a = R' * a;    % 3 x n
    f = R' * f;    % 3 x n
    f = f + X(1,:)';   % shift; implicit broadcasting of X' (3x1) over columns

    % Output: [radius; a(1:3,:); f(1:3,:)]  -> 7 x n
    cyl = [radius; a ./ norma; f];
end

function dist = cyldist(cyl, X, pos)
    % cyldist(cyl,X,pos) computes distances of points from cylinder
    if nargin < 3, pos = 1; end
    r = cyl(1);
    a = cyl(2:4);
    f = cyl(5:7);

    xf = X - f';
    nxf = sum(xf.^2, 2);
    if pos == 1
        dist = abs(sqrt(nxf - (xf*a).^2) - r);
    else
        dist = sqrt(nxf - (xf*a).^2) - r;
    end
end

function [s, r, c] = acfit(X, a, sonly)
    if nargin < 3, sonly = false; end
    n = size(X,1);
    t = a(1); l = a(2);

    % Projection direction and plane basis
    au = [sin(t)*cos(l), sin(t)*sin(l), cos(t)];
    ath = [cos(t)*cos(l), cos(t)*sin(l), -sin(t)];
    ala = [-sin(l), cos(l), 0];

    % Projected coordinates
    x = X * ala';
    y = X * ath';
    Xs = [x, y];

    % Fit circle
    M = [-2*Xs, ones(n,1)];
    h = -sum(Xs.^2, 2);
    p = M \ h;
    r = sqrt(p(1)^2 + p(2)^2 - p(3));
    cs = p(1:2);
    TR = [ala; ath; au]';
    c = TR * [cs; 0];
    s = sum((M*p - h).^2);

    if sonly
        s = s;
    end
end

function [c, r, a] = cylinder(X)
    % cylinder(X): fits a right circular cylinder into the points X
    Xm = mean(X, 1);
    Xs = X - Xm;
    n = size(X, 1);

    % Cost function
    sqsum = @(x) acfit(X, x, true);
    x0 = [0.1, 0.1];
    % opts = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton');
    % [x, ~] = fminunc(sqsum, x0, opts);
    [x, ~] = fminunc(sqsum, x0);

    a = x;
    au = [sin(x(1))*cos(x(2)), sin(x(1))*sin(x(2)), cos(x(1))];
    fprintf('axis unit vector: [ %.5f, %.5f, %.5f ]\n', au(1), au(2), au(3));

    ath = [cos(x(1))*cos(x(2)), cos(x(1))*sin(x(2)), -sin(x(1))];
    ala = [-sin(x(2)), cos(x(2)), 0];

    % Projected coordinates
    xproj = X * ala';
    yproj = X * ath';

    [~, r, c] = acfit(Xs, a);
    c = c' + Xm;
end


