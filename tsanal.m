%% time series analysis of structural monitoring data
% Field course of Structural Geodesy, measurement at 2018.10.11

clear all; close all

function sec = hms2s(hms)  % seconds of the day
    sp = strsplit(hms,':');
    h = str2num(sp{1});
    m = str2num(sp{2});
    s = str2num(sp{3});
    sec = 3600*h+60*m+s;
end

% read GNSS data
% P.No  Y   X    h    code    date   time
[PNo, Y, X, h, code, dat, time1] = textread('181011G.txt', '%s %f %f %f %s %s %s', 'headerlines',7);

% seconds from the first epoch
n1 = length(h);
s0 = hms2s(time1{1});
t1 = zeros(1,n1);
for i=2:n1
    t1(i) = hms2s(time1{i})-s0;
end

% are there missing epochs?
dt1 = diff(t1);
me = length(dt1(dt1>1));
if me > 0
    fprintf('there are missing epochs\n')
    [mx,ix] = max(dt1);
    fprintf('max: %d index: %d\n',mx,ix)
    % linear interpolation of missing data
    tmax1 = t1(end); 
    length(t1); 
    tf1 = 1:tmax1;
    hi = interp1(t1,h,tf1);
    % save interpolated values
    f1 = fopen('G.dat','w');
    for i=1:tmax1
        fprintf(f1,"%d %.3f\n",tf1(i),hi(i));
    end
    fclose(f1);
else
    hi = h;
    tf1 = t1;
end


% plot heights
figure(1);
plot(t1,h)
xlim([0, 1700])
xlabel('time (s)');
ylabel('height (m)');
title('GPS bridge monitoring 2018.10.11.');

% calculate autocovariance function
pkg load signal  % Octave
[Rh, lag] = xcorr(hi-mean(hi),'unbiased');
% plot
figure(2);
plot(lag,1e4*Rh)  % convert to cm^2
xlabel('time (s)');
ylabel('cm^2')
xlim([0, 1700])
title('GPS bridge height ACF');
hold on
plot([0,1700],[0,0],'r-')

%% determine amplitude and phase of a T period signal component
% subtract sine function with shift, amplitude and calculate SSR
function ssr=lossnv(dt,A,T,data)  % non-vectorized version
    nd = length(data);
    t = 0:(nd-1);
    ys = A*sin(2*pi/(T)*(t+dt));
    % ssr = log10(sum((data-ys).^2)/sum(data.^2));
    ssr = sum((data-ys).^2)/sum(data.^2);
end

function ssr = loss(dt, A, T, data)  % vectorized version
    nd = length(data);
    t = 0:(nd-1);

    % Normalization factor
    denom = sum(data.^2);

    % Reshape for broadcasting
    dt3   = reshape(dt, [], 1, 1);     % Nd × 1 × 1
    A3    = reshape(A, 1, [], 1);      % 1 × Na × 1
    t3    = reshape(t, 1, 1, []);      % 1 × 1 × Nt
    data3 = reshape(data, 1, 1, []);   % 1 × 1 × Nt

    % Predicted signals (all dt × all A × time)
    ys = A3 .* sin(2*pi/T .* (t3 + dt3));

    % SSR for each (dt, A)
    ssr = squeeze(sum((data3 - ys).^2, 3)) / denom;
end

% find and remove T=680s sine wave
nd = length(hi);
T = 680.0;
dt = 0:(T/100):(T);
A = 0.002:0.0002:0.02;
t = 0:(nd-1);
hic = hi-mean(hi);
ssr = loss(dt,A,T,hic);
[minc,i] = min(ssr);
[totalmin,j] = min(minc);
totalmax = max(max(ssr));
Amin = A(j);
dtmin = dt(i(j));
fprintf('best sinusoid %.1f sec: A: %.3f dt: %.1f, min. ssr: %.5f\n', ...
        T, Amin,dtmin,totalmin)
% best sinusoid 680.0 sec: A: 0.009 dt: 476.0, min. ssr: 0.83949
figure(3)
contour(A,dt,ssr,(totalmin+0.01):0.05:totalmax)
ylabel('dt (s)')
xlabel('A (m)')

% subtract 680 sec best sinusoid
hic = hic - 0.009*sin(2*pi/(680)*(t'+476));

% plot ACF
[Rh, lag] = xcorr(hic-mean(hic),'unbiased');
% plot
figure(4);
plot(lag,1e4*Rh)  % convert to cm^2
xlabel('time (s)');
ylabel('cm^2')
xlim([0, 1700])
title('GPS bridge height ACF after removing 680 s period sine');
hold on
plot([0,1700],[0,0],'r-')

% find and remove T=90 s sine wave
T = 90;
dt = 0:(T/100):(T);
A = 0.002:0.0002:0.02;
t = 0:(nd-1);
hic = hic-mean(hic);
ssr = loss(dt,A,T,hic);
[minc,i] = min(ssr);
[totalmin,j] = min(minc);
totalmax = max(max(ssr));
Amin = A(j);
dtmin = dt(i(j));
fprintf('best sinusoid %.1f sec: A: %.3f dt: %.1f, min. ssr: %.5f\n', ...
        T, Amin,dtmin,totalmin)
% best sinusoid 90.0 sec: A: 0.010 dt: 16.2, min. ssr: 0.71307
figure(5)
contour(A,dt,ssr,(totalmin+0.01):0.05:totalmax)
ylabel('dt (s)')
xlabel('A (m)')

% subtract 90 sec best sinusoid
hic = hic - 0.010*sin(2*pi/(90)*(t'+16.2));

% plot ACF
[Rh, lag] = xcorr(hic-mean(hic),'unbiased');
% plot
figure(6);
plot(lag,1e4*Rh)  % convert to cm^2
xlabel('time (s)');
ylabel('cm^2')
xlim([0, 1700])
title('GPS bridge height ACF after removing 90 s period sine');
hold on
plot([0,1700],[0,0],'r-')

% plot filtered signal
figure(1)
hold on
plot(t,hic+mean(hi),'r-')
xlim([0, 1700])
grid on

% percent std decrease
pct = 100*(std(hi)-std(hic+mean(hi)))/std(hi);  % 22.635

% read and plot bus crossing times
[time] = textread('bus.txt', '%s');
% seconds from the first epoch
n = length(time);
tb = zeros(1,n);
for i=1:n
    tb(i) = hms2s(time{i})-s0;
end
nb = length(tb);
i1 = 1:12;  % articulated bus
i2 = 13:19;  % tour bus
i3 = 20:nb;  % simple bus

figure(1)
scatter(tb(i1),113.87*ones(length(tb(i1)),1),25,'k','filled')
scatter(tb(i2),113.87*ones(length(tb(i2)),1),25,'m','filled')
scatter(tb(i3),113.87*ones(length(tb(i3)),1),25,'b','filled')
grid on
legend('GPS','filtered','articulated bus','tour bus','simple bus')





