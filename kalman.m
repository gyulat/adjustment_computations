% Kalman filter step
% x0, P0: input state and covariance
% Bu    : control
% z, R  : measurement and its noise covariance
% Q     : state noise covariance
% F     : state transition matrix
% H     : measurement matrix
% x1, P1: updated state and covariance

function [x1,P1] = kalman(x0,P0,Bu,z,Q,R,F,H)
      
      xpr = F*x0 + Bu;  % prediction
      Ppr = F*P0*F' + Q;  % prediction covariance
      if (isempty(z)), % predict only, 
          xup = xpr;
          x1 = xpr; P1 = Ppr;
          return
      end;
      S=H*Ppr*H'+R; % covariance of innovation
      K=Ppr*H'/S;   % Kalman gain
      y=z-H*xpr;    % innovation
      xup = xpr;
      x1=xpr+K*y;  % measurement update
      P1=Ppr-K*H*Ppr; 
end
