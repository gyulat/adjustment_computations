%% autocovariance function of sines
clear all; close all
% signal samples
ns = 1200;
% create sine functions with shift and periodicity
Tp1 = 700;
a1 = 0.5;
ts1 = 50;
Tp2 = 90;
a2 = 1;
ts2 = 0;
t = 0:(ns-1);
ys1 = a1*sin(2*pi/(Tp1)*(t+ts1));
ys2 = a2*sin(2*pi/(Tp2)*(t+ts2));

% sum of sines
ys = ys1 + ys2;

noise = 1;
snoise ='';

if noise
    % white noise
    yn = 0.5*randn(1,ns);
    % red noise
    alpha = 0.8;
    yr = zeros(1,ns);
    yr(1) = yn(1);
    for i=2:ns
        yr(i) = alpha*yr(i-1)+yn(i);
    end
    ys = ys + yr;
    snoise = 'with red noise';
end

pkg load signal % Octave

figure(1)
plot(t,ys)
hold on
plot(t,ys1)
plot(t,ys2)
grid on
[ryy,lag]=xcorr(ys-mean(ys),'unbiased');
figure(2)
plot(lag,ryy);
axis([0,1200])
grid on
title(sprintf('sum of sines autocovariance %s', snoise))
xlabel('lag (s)')


if noise
    % plot red noise ACF
    [ryyr,lagr]=xcorr(yr-mean(yr),'unbiased');
    figure(3)
    plot(lagr,ryyr);
    axis([0,1200])
    grid on
    title('red noise autocovariance')
    xlabel('lag (s)')
end


