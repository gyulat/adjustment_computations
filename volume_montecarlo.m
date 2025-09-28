% volumetric data error calculation with Monte-Carlo
clc; clear all; close all
% load measurements

data = load('volume.txt');
x = data(:,1);
y = data(:,2);
z1 = data(:,3);
z2 = data(:,4);
stdev = data(:,5);
ndata = length(x);

% vectorized volume calculation for both epochs
vol1 = volume(x,y,z1);
vol2 = volume(x,y,z2);
vol12 = vol1-vol2;
fprintf('calculated volume: %.2f m³\n',vol12)


% Monte Carlo calculation of the volume error
Nsim = 10000;
volerr = zeros(Nsim,1);
dat1 = [x,y,zeros(ndata,1)];
dat2 = [x,y,zeros(ndata,1)];
for j=1:Nsim
    zerr1 = z1+stdev.*randn(ndata,1);
    zerr2 = z2+stdev.*randn(ndata,1);
    vol1 = volume(x,y,zerr1);
    vol2 = volume(x,y,zerr2);
    volerr(j) = vol1-vol2 - vol12;
end

printf('volume error std from %d runs: %.2f m³\n', Nsim, std(volerr))
#    calculated volume: 35288.88 m³
#    volume error std from 10000 runs: 18.77 m³

figure(1)
hist(volerr, 30)
title(sprintf('Volume error distribution for %d runs',Nsim))
xlabel('m³')
