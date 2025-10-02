function kalman_demo
    % Simulate
    duration = 100;
    dt = 0.5;
    [pos, posmeas, poshat] = kalman(duration, dt);

    % Time vector
    t = 0:dt:(duration-dt);

    % Plot the results
    figure;
    plot(t,posmeas,'DisplayName','measured'); hold on;
    plot(t,pos,'DisplayName','true');
    plot(t,poshat,'DisplayName','estimated');
    legend('Location','best');
    xlabel('time (s)');
    ylabel('position (m)');
    title('Kalman filter simulation');
    grid on;
    hold off;
end

function [pos,posmeas,poshat] = kalman(duration,dt,accelnoise)
    % Kalman filter simulation
    % duration = simulation length (s)
    % dt = timestep (s)
    % accelnoise = acceleration noise (m/s^2)

    if nargin < 3
        accelnoise = 0.2; % default
    end

    measnoise  = 3;   % position measurement noise (m)

    % State transition and measurement matrices
    F = [1 dt; 0 1];
    H = [1 0];

    % Initial state
    x = [0;0];
    xhat = x;

    % Process noise covariance
    Q = accelnoise^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2];

    % Initial estimation covariance
    P = Q;

    % Measurement covariance
    R = measnoise^2;

    % Number of epochs
    dimt = floor(duration/dt);

    % Storage
    pos     = zeros(1,dimt);
    posmeas = zeros(1,dimt);
    poshat  = zeros(1,dimt);

    % Main loop
    for k = 1:dimt
        % Simulate process
        ProcessNoise = accelnoise * [(dt^2/2)*randn; dt*randn];
        x = F*x + ProcessNoise;

        % Simulate measurement
        MeasNoise = measnoise * randn;
        z = H*x + MeasNoise;

        % Innovation
        Inn = z - H*xhat;

        % Innovation covariance
        S = H*P*H' + R;

        % Kalman gain
        K = (F*P*H') / S;

        % State estimate
        xhat = F*xhat + K*Inn;

        % Covariance update
        P = F*P*F' + Q - (K*H*P*F');

        % Save results
        pos(k)     = x(1);
        posmeas(k) = z;
        poshat(k)  = xhat(1);
    end
end

