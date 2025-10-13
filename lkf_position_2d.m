% function for Kalman filtering to estimate the trajectory of an object moving in a plane
function [xy_new, P] = lkf_position_2d(t,xy,q,r)
    % t - (N,1) vector of time epochs
    % xy - (N,2) matrix of x,y measurements
    % q, r - (4,1) and (2,1) vectors of diagonal elements of Q and R matrices
    % initialize filter
    x0 = [xy(1,1);xy(1,2);0;0]; % zero initial velocity 
    Q = diag(q);
    P0 = 10*Q;
    R = diag(r);
    N = length(t);
    H = eye(2,4);
    Bu = zeros(4,1);
    xy_new = zeros(N,2);
    xy_new(1,:) = x0(1:2);
    P = zeros(N,16);
    P(1,:) = P0(:)';
    % estimate x,y
    for i=2:N
        dt = t(i)-t(i-1);
        F = eye(4); F(1,3) = dt; F(2,4) = dt;
        z = xy(i,:)';
        [x1,Pi] = kalman(x0,P0,Bu,z,Q,R,F,H);
        xy_new(i,:) = x1(1:2)';
        x0 = x1; % update state
        P(i,:) = Pi(:)';
    end
end 
    
function [x1,P1] = kalman(x0,P0,Bu,z,Q,R,F,H)    
      xpr = F*x0 + Bu;  % prediction
      Ppr = F*P0*F' + Q;  % prediction covariance
      S=H*Ppr*H'+R; % covariance of innovation
      K=Ppr*H'/S;   % Kalman gain
      y=z-H*xpr;    % innovation
      x1=xpr+K*y;   % measurement update
      P1=Ppr-K*H*Ppr; 
end
