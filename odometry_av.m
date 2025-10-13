function main
    clc; clear all; close all
    % Load data
    data = load('a.txt');
    t  = data(:,1);
    a1 = data(:,2);
    a2 = data(:,3);
    av = load('om_interp.txt');
    om = av(:,2);  % same time epochs

    % Plot acceleration and angular velocity data
    figure;
    plot(t,a1,'r-','LineWidth',1,'DisplayName','acceleration 1, m/s^2'); hold on;
    plot(t,a2,'b-','LineWidth',1,'DisplayName','acceleration 2, m/s^2');
    plot(t,om,'m-','LineWidth',1,'DisplayName','angular velocity, rad/s');
    grid on;
    title('Galaxy S2 sensor data','FontWeight','bold');
    xlabel('time (s)');
    ylabel('acceleration (m/s^2) / angular velocity (rad/s)');
    legend('Location','northeast');
    hold off;


    % Run EKF wheel odometry
    xe = EKFwheelav(t,a1,a2,om,0.17,5.0,0.5);

    % Plot EKF results
    figure;
    plot(xe(:,1),xe(:,2),'b-','DisplayName','EKF distance, m'); hold on;
    plot(xe(:,1),xe(:,3),'g-','DisplayName','EKF velocity, m/s');
    plot(xe(:,1),xe(:,4),'r-','DisplayName','EKF acceleration, m/s^2');
    plot(xe(:,1),xe(:,5),'Color',[1 0.5 0],'DisplayName','angular position, rad');
    %plot(xe(:,1),xe(:,6),'b--','DisplayName','model distance, m');
    %plot(xe(:,1),xe(:,7),'g--','DisplayName','model velocity, m/s');
    %plot(xe(:,1),xe(:,8),'r--','DisplayName','model acceleration, m/s^2');
    xlim([1,11]);
    grid on;
    title('Results of the Extended Kalman filter (EKF)','FontWeight','bold');
    xlabel('time (s)');
    legend('Location','northwest');
    hold off;
    
    % compare with distances without gyro measurements
    dod = load('dist_odometry.txt');
    dist = dod(:,2);
    dd = xe(:,2)-dist;
    figure;
    plot(xe(:,1),dd,'b-','DisplayName','EKF distance differences, m'); hold on;
    xlim([1,11]);
    grid on;
    title('Distance differences of EKF with gyro - without gyro','FontWeight','bold');
    xlabel('time (s)');
    legend('Location','northwest');
    hold off;
    
end

%% State propagation function
function xnext = f(x,dt)
    f1 = x(1) + x(2)*dt + 0.5*x(3)*dt^2;
    f2 = x(2) + x(3)*dt;
    f3 = x(3);
    xnext = [f1; f2; f3];
end

%% Measurement function
function z = h(x,rs,rw)
    g = 9.81;
    h1 = -g*sin(x(1)/rw) + x(3)*cos(x(1)/rw) - x(3)*rs/rw;
    h2 = -g*cos(x(1)/rw) - x(3)*sin(x(1)/rw) - (x(2))^2*rs/(rw^2);
    h3 = -x(2)/rw;
    z = [h1; h2; h3];
end

%% Jacobian of measurement function
function Hk = Hjac(x,rs,rw)
    g = 9.81;
    Hk = zeros(3,3);
    Hk(1,1) = -g/rw*cos(x(1)/rw) - x(3)/rw*sin(x(1)/rw);
    Hk(1,2) = 0.0;
    Hk(1,3) = cos(x(1)/rw) - rs/rw;
    Hk(2,1) = g/rw*sin(x(1)/rw) - x(3)/rw*cos(x(1)/rw);
    Hk(2,2) = -2*x(2)*rs/(rw^2);
    Hk(2,3) = -sin(x(1)/rw);
    Hk(3,1) = 0.0;
    Hk(3,2) = -1.0/rw;
    Hk(3,3) = 0.0;
end

%% One EKF step
function [xkm,xk1,Pk1] = EKFstep(xk,Pk,A,H,Q,R,zk,dt,rs,rw)
    % Prediction
    xkm = f(xk,dt);
    Pkm = A*Pk*A' + Q;

    % Kalman gain
    Kk = Pkm*H' / (H*Pkm*H' + R);

    % Update
    xk1 = xkm + Kk*(zk - h(xkm,rs,rw));
    Pk1 = Pkm - Kk*H*Pkm;
end

%% EKF wheel odometry
function xe = EKFwheelav(t,a1,a2,om,q,r,rom)
    rs = 0.095;   % sensor offset (m)
    rw = 0.35;    % wheel radius (m)
    nt = length(t);

    % xe = zeros(nt,5);
    xe = zeros(nt,8);
    xe(:,1) = t;

    % Initialization
    Q = q^2*eye(3);
    R = r^2*eye(3); R(3,3) = rom^2; % angular velocity noise variance
    xk = zeros(3,1);
    Pk = q^2*diag([0,0,1]);

    for i = 2:nt
        dt = t(i)-t(i-1);
        Hk = Hjac(xk,rs,rw);
        A  = eye(3) + [0 dt 0.5*dt^2; 0 0 dt; 0 0 0];
        zk = [a1(i); a2(i); om(i)];

        [xkm1,xk1,Pk1] = EKFstep(xk,Pk,A,Hk,Q,R,zk,dt,rs,rw);
        xe(i,2:4) = xk1';
        xe(i,6:8) = xkm1;
        xk = xk1;
        Pk = Pk1;

        ome = xk1(1)/rw; % angular position
        xe(i,5) = atan2(sin(ome),cos(ome)); % wrap to [-pi,pi]
    end
end




