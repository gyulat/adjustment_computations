%% Fitting cicles to chimney sections
clear all; close all

dat=load('chimney_sections.txt');
x=dat(:,1);
y=dat(:,2);
z=dat(:,3);
%figure(1)
%scatter3(x,y,z,10,'bo',"filled")
%axis equal
nd = length(z);
fprintf('number of datapoints: %d\n',nd)
tol = 0.5; % height tolerance to section selection


%% determine sections
dz = diff(z);
ns = sum(abs(dz)>tol) + 1;  % max number of sections
maxs = 10;  % initial maximum number of points in a section
ix = zeros(ns,maxs); % point indexes of sections
hs = zeros(ns,1); % heights of sections
nk = zeros(ns,1); % number of points within a section
j=1;  % number of section
k=1;  % number of point inside section
ix(j,1) = 1; % first point
hs(j) = z(j);
nk(1) = 1;
sfound = j; % total number of sections found
oldsection = false;
for i=2:nd
    oldsection = false;
    % check existing section heights
    for m = 1:sfound %-1
        if abs(z(i)-hs(m))<tol  % add point to existing section
            oldsection = true;
            j = m;
            k = nk(m) + 1;
            ix(m,k) = i; % next point in section
            nk(m) = k;
            break
        end  
    end
    % not in previous sections 
    if ~oldsection
        % oldsection = false;
        if abs(z(i)-hs(j))<tol  % add another point to new section
            k = k+1;
            ix(j,k) = i;
            nk(j) = k;
        else
            % start new section
            j = sfound + 1;
            sfound = j; % number of sections found
            nk(j) = 1;
            k = 1;
            ix(j,k) = i; % first point in section
            hs(j) = z(i); % set new section height
        end
    end
end
% chop zero parts
ix = ix(1:sfound,:);
nk = nk(1:sfound);
hs = hs(1:sfound);
fprintf('number of sections found: %d\n',sfound)
fprintf('number of datapoints in all sections: %d\n',sum(nk))

% plot sections
figure(1);
hold on
for i = 1:sfound
    % select i-th section
    xi = x(ix(i,1:nk(i)));
    yi = y(ix(i,1:nk(i)));
    zi = z(ix(i,1:nk(i)));
    scatter3(xi,yi,zi,10,'o',"filled")
end
axis equal

% fit circles to each chimney section
par = zeros(sfound,3); % fitted circle parameters
for i = 1:sfound
    % select i-th section
    xi = x(ix(i,1:nk(i)));
    yi = y(ix(i,1:nk(i)));
    % LSQ circle fit
    A = [xi,yi,ones(length(xi),1)]; b = -[xi.^2+yi.^2]; p = A\b;
    xc = -0.5*p(1); yc = -0.5*p(2);
    R = sqrt((p(1)^2+p(2)^2)/4-p(3));
    par(i,:) = [xc,yc,R];
end

% plot circle centers and radii
figure(2)
plot(par(:,1),par(:,2),'b-')
hold on
plot(par(1,1),par(1,2),'rd')
axis equal
grid on
title('chimney axis points')
figure(3)
plot(hs,par(:,3),'r-')
grid on
xlabel('height (m)')
ylabel('radius (m)')
title('chimney section radii')
