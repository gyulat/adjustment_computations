% To run this script, you need the following data available in the MATLAB
% workspace:
% * pos_xy - original, real position data
% * pos_xy_m - 'measured' position with noise
% * q, r - arbitrary vectors
clc; clear all; close all;

data = load('trajectory.txt');
pos_t = data(:,1);
pos_xy_m = data(:,2:3);
pos_x = data(:,4); pos_y = data(:,5);

% Run LKF 
% system state noise variance
q = 0.01*ones(1,4);   % position
q(3:4) = 0.01*[1,1];  % velocity
r = 0.1*ones(1,2);
[pos_xy, P] = lkf_position_2d(pos_t, pos_xy_m, q, r);

% Extract
fpos_x = pos_xy(:,1); fpos_y = pos_xy(:,2);
rpos_x = pos_xy_m(:,1); rpos_y = pos_xy_m(:,2);

% Compute MSE for residuals
get_mse = @(x,y) 1/length(x) * sum((x-y).^2);
mse_x_ns = get_mse(pos_x, rpos_x);
mse_x_kf = get_mse(pos_x, fpos_x);
mse_y_ns = get_mse(pos_y, rpos_y);
mse_y_kf = get_mse(pos_y, fpos_y);

fprintf('q coordinates: %.3f m^2/s^2\n',q(1));
fprintf('q velocities: %.3f m^4/s^4\n',q(3));
fprintf('MSE xmeas: %.4f ymeas: %.4f\n', mse_x_ns, mse_y_ns);
fprintf('MSE xfilt: %.4f yfilt: %.4f\n', mse_x_kf, mse_y_kf);

% Display the results
figure(1, 'name', 'Kalman filtering test results: coordinates vs. time');

subplot(211); plot(pos_t, rpos_x, '--r'); hold on;
plot(pos_t, fpos_x, 'b'); plot(pos_t, pos_x, 'k');
legend(['Measured, MSE=' num2str(mse_x_ns)], ...
    ['Filtered, MSE=' num2str(mse_x_kf)], 'Real');
ylabel('x');

subplot(212); plot(pos_t, rpos_y, '--r'); hold on;
plot(pos_t, fpos_y, 'b'); plot(pos_t, pos_y, 'k');
legend(['Measured, MSE=' num2str(mse_y_ns)], ...
    ['Filtered, MSE=' num2str(mse_y_kf)], 'Real');
ylabel('y');
xlabel('t');

% Plot tracks
figure(2, 'name', 'Kalman filtering test results: trajectories');
plot(pos_xy_m(:,1), pos_xy_m(:,2),'-*r');
hold on;
plot(fpos_x, fpos_y,'-ob');
plot(pos_x, pos_y,'k');
legend('Measured', 'Filtered', 'Ideal');
