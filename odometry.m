function main
    % Load data
    data = load('a.txt');
    t  = data(:,1);
    a1 = data(:,2);
    a2 = data(:,3);

    % Plot acceleration data
    figure;
    plot(t,a1,'r-','LineWidth',1,'DisplayName','acceleration 1, m/s^2'); hold on;
    plot(t,a2,'b-','LineWidth',1,'DisplayName','acceleration 2, m/s^2');
    grid on;
    title('Galaxy S2 sensor data','FontWeight','bold');
    xlabel('time (s)');
    ylabel('acceleration (m/s^2)');
    legend('Location','northeast');
    hold off;

    % Run EKF wheel odometry
    xe = EKFwheel(t,a1,a2,0.17,5.0);

    % Plot EKF results
    figure;
    plot(xe(:,1),xe(:,2),'b-','DisplayName','EKF distance, m'); hold on;
    plot(xe(:,1),xe(:,3),'g-','DisplayName','EKF velocity, m/s');
    plot(xe(:,1),xe(:,4),'r-','DisplayName','EKF acceleration, m/s^2');
    plot(xe(:,1),xe(:,5),'Color',[1 0.5 0],'DisplayName','angular position, rad');
    xlim([1,11]);
    grid on;
    title('Results of the Extended Kalman filter (EKF)','FontWeight','bold');
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
    z = [h1; h2];
end

%% Jacobian of measurement function
function Hk = Hjac(x,rs,rw)
    g = 9.81;
    Hk = zeros(2,3);
    Hk(1,1) = -g/rw*cos(x(1)/rw) - x(3)/rw*sin(x(1)/rw);
    Hk(1,2) = 0.0;
    Hk(1,3) = cos(x(1)/rw) - rs/rw;
    Hk(2,1) = g/rw*sin(x(1)/rw) - x(3)/rw*cos(x(1)/rw);
    Hk(2,2) = -2*x(2)*rs/(rw^2);
    Hk(2,3) = -sin(x(1)/rw);
end

%% One EKF step
function [xk1,Pk1] = EKFstep(xk,Pk,A,H,Q,R,zk,dt,rs,rw)
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
function xe = EKFwheel(t,a1,a2,q,r)
    rs = 0.095;   % sensor offset (m)
    rw = 0.35;    % wheel radius (m)
    nt = length(t);

    xe = zeros(nt,5);
    xe(:,1) = t;

    % Initialization
    Q = q^2*eye(3);
    R = r^2*eye(2);
    xk = zeros(3,1);
    Pk = q^2*diag([0,0,1]);

    for i = 2:nt
        dt = t(i)-t(i-1);
        Hk = Hjac(xk,rs,rw);
        A  = eye(3) + [0 dt 0.5*dt^2; 0 0 dt; 0 0 0];
        zk = [a1(i); a2(i)];

        [xk1,Pk1] = EKFstep(xk,Pk,A,Hk,Q,R,zk,dt,rs,rw);
        xe(i,2:4) = xk1';
        xk = xk1;
        Pk = Pk1;

        om = xk1(1)/rw; % angular position
        xe(i,5) = atan2(sin(om),cos(om)); % wrap to [-pi,pi]
    end
end




