clear all; close all

pkg load statistics % Octave
pconf = 0.95;
fprintf('\n*** Montsalvens deformation detection network ***\n\n')

% first epoch results
load montsalvens1.mat
m0_apost1 = m_0_aposteriori;
ssq1 = sum_of_squares;
f1 = observations-unknowns+network_defect;
x1 = x; A1 = A; b1 = b; P1 = P;
N1 = N; n1 = n;
Cll1 = C_ll;
Qxx1 = inv(N1); % chol(Qxx1); % not positive definite
Qxx1 = Qxx1(1:unknowns, 1:unknowns); % (29,29)
QUU1 = A1*Qxx1*A1';
Qvv1 = Cll1 - QUU1;
v1 = full(A1*x1-b1);
R1 = eye(observations)-QUU1*P1;
XYZ1 = XYZ;
Cxx1 = C_xx;

% second epoch results
load montsalvens2.mat
m0_apost2 = m_0_aposteriori;
ssq2 = sum_of_squares;
f2 = observations-unknowns+network_defect;
x2 = x; A2 = A; b2 = b; P2 = P;
N2 = N; n2 = n;
Cll2 = C_ll;
Qxx2 = inv(N2);
Qxx2 = Qxx2(1:unknowns, 1:unknowns);
QUU2 = A2*Qxx2*A2';
Qvv2 = Cll2 - QUU2;
v2 = full(A2*x2-b2);
R2 = eye(observations)-QUU2*P2;
XYZ2 = XYZ;
Cxx2 = C_xx;

% save workspace for plotting
X1 = XYZ1(:,1); X2 = XYZ2(:,1); 
Y1 = XYZ1(:,2); Y2 = XYZ2(:,2); 
Cov1 = Cxx1; Cov2 = Cxx2;
edg = [1,2;1,5;1,6;1,8;1,9;1,10;2,6;2,8;2,5;2,4;2,9;2,10;2,11;...
        3,2;3,5;3,6;3,10;3,11;3,12;3,7;4,3;4,1;4,5;4,6;4,9;4,11;4,12;5,6]; 
stable = [1:7];

save -mat deform2d.mat
 
% return;

% between epochs 1 and 2
dX = XYZ2(:,1)-XYZ1(:,1);
dY = XYZ2(:,2)-XYZ1(:,2);
fprintf('Epochs 2-1 \n')
fprintf('Point  dX(mm) dY(mm) \n')
for i=1:length(XYZ1)
    fprintf('  %2s  %6.2f %6.2f\n',Points{i},1000*dX(i),1000*dY(i));
end

% Congruency test for all the reference points
% ref. points 1,2,3,4,6,7,9
ptol = 1.0e-12; % tolerance for pseudoinverses - very important
rpi = 1:7; % ref. point indexes in Points
m = length(rpi); % number of reference points
icr = Indexes(rpi,1:2);  % x y coord. indexes
icrf = icr(:); % flatten
icq = [1:2:2*m,2:2:2*m];
Qd = Qxx1(icq,icq)+Qxx2(icq,icq);
D = x1(icrf)-x2(icrf);
qd = D'*pinv(Qd,ptol)*D  % 54.759
fd = length(icrf) - network_defect;  % 11

q = ssq1 + ssq2;  % 61.989 - same v1'*P1*v1 + v2'*P2*v2
f = f1 + f2;
s02 = q/f;  % 1.0688
T = qd/fd / (q/f)  %  4.6578
% Fisher F distribution inverse CDF
Tcrit = finv(pconf,fd,f)  %  1.9580
if T > Tcrit
    fprintf('*** not all reference points are stable at level %.2f\n',pconf)
else
    fprintf('+++ the stability of all reference points cannot be rejected at level %.2f\n',pconf)
end

% pointwise test for reference point stability 
% 
% reference points
Pdc = pinv(Qd,ptol);
Dc = D;
maxq = 0;
for i=1:m % m
    % split of Dc, Qdc for point i
    % for reference points
    Ind = reshape(1:2*m,m,2);
    ip = Ind(i,:);
    noti = (1:m)(setdiff(1:end,i));
    in = Ind(noti,:);
    ifp = ip(:);
    ifn = in(:);
    Dp = Dc(ifp);
    Dn = Dc(ifn);
    Pnn = Pdc(ifn,ifn); 
    Ppn = Pdc(ifp,ifn); 
    Pnp = Pdc(ifn,ifp);
    Ppp = Pdc(ifp,ifp);
    Pdf = [Pnn,Pnp;Ppn,Ppp];
    Dpo = pinv(Ppp)*Ppn*Dn + Dp;
    qdp = Dpo'*Ppp*Dpo;
    Qnninv = Pnn-Pnp*pinv(Ppp,ptol)*Ppn;
    qdn = Dn'*Qnninv*Dn;
    if qdp > maxq
        maxq = qdp;
        maxi = i;
    end
    fprintf('   q_%s: %.3f\n',Points{i},qdp)   
end
fprintf('max. q is at point %s: %.3f\n',Points{i}, maxq)  % 

% New congruency test without reference point 3.
qd1 = qd - maxq;  % 10.634
fd1 = fd - 2;  % 9
T = qd1/fd1 / (q/f);  %  1.1055
% Fisher F distribution inverse CDF
Tcrit = finv(pconf,fd1,f);  %  2.0458
if T > Tcrit
    fprintf('*** not all reference points are stable at level %.2f\n',pconf)
else
    fprintf('+++ the stability of all reference points cannot be rejected at level %.2f\n',pconf)
end

%    *** not all reference points are stable at level 0.95
%       q_1: 1.297
%       q_2: 3.810
%       q_3: 44.109
%       q_4: 16.281
%       q_6: 1.059
%       q_7: 2.228
%       q_9: 0.018
%    max. q is at point 9: 44.109
%    +++ the stability of all reference points cannot be rejected at level 0.95


%%  S-transformation to the datum of stable reference points
npts = length(XYZ1);  % number of points

% constrain coordinates of the 6 stable reference points 
% indexes of fixed coordinates
istable = [1,2,4,5,6,7];
d = zeros(1,2*npts);
ip = fix(Indexes(:,1)/2);
for i=1:npts
    if any (istable == ip(i))
        d(i) = 1;  % for x
        d(i+12) = 1;  % for y
    end
end
% datum selector matrix E
E = diag(d); % (24,24)

% set up G-matrix for epoch 1
G = zeros(2*npts,3); % (24,3)
% separate coordinates
X = XYZ1(:,1); Y = XYZ1(:,2);
% indexes of points
ixsol = reshape(Indexes(:,1:2),2*npts,1);
ix = Indexes(:,1); iy = Indexes(:,2);
ip = fix(ix/2);
G(1:12,1) = 1; G(13:end,2) = 1;
G(1:12,3) = Y(ip); G(13:end,3) = -X(ip);
% S-matrix
S = eye(2*npts) - G*inv(G'*E*G)*G'*E;
% solution
xS1 = S*x1(ixsol);
% adjusted coordinates
XS1 = XYZ_0(:,1)+0.001*xS1(1:12);
YS1 = XYZ_0(:,2)+0.001*xS1(13:end);
if false
fprintf('Epoch 1.\n')
fprintf('Point      X          Y \n')
for i=1:npts
    fprintf('  %2s  %10.4f %10.4f\n',Points{i},XS1(i),YS1(i));
end
end
% covariance and cofactor mx.
CS1 = S*Cxx1(1:2*npts,1:2*npts)*S';
QS1 = S*Qxx1(1:2*npts,1:2*npts)*S';

% set up G-matrix for epoch 2
G = zeros(2*npts,3); % (24,3)
% separate coordinates
X = XYZ2(:,1); Y = XYZ2(:,2);
% indexes of points
ixsol = reshape(Indexes(:,1:2),2*npts,1);
ix = Indexes(:,1); iy = Indexes(:,2);
ip = fix(ix/2);
G(1:12,1) = 1; G(13:end,2) = 1;
G(1:12,3) = Y(ip); G(13:end,3) = -X(ip);
% S-matrix
S = eye(2*npts) - G*inv(G'*E*G)*G'*E;
% solution
xS2 = S*x2(ixsol);
% adjusted coordinates
XS2 = XYZ_0(:,1)+0.001*xS2(1:12);
YS2 = XYZ_0(:,2)+0.001*xS2(13:end);

if false
fprintf('Epoch 2.\n')
fprintf('Point      X          Y \n')
for i=1:npts
    fprintf('  %2s  %10.4f %10.4f\n',Points{i},XS2(i),YS2(i));
end
end
% covariance and cofactor mx.
CS2 = S*Cxx2(1:2*npts,1:2*npts)*S';
QS2 = S*Qxx2(1:2*npts,1:2*npts)*S';

% save workspace for plotting
X1 = XS1; X2 = XS2; 
Y1 = YS1; Y2 = YS2; 
Cov1 = CS1; Cov2 = CS2;
edg = []; 
stable = istable;

save -mat deform2d_6ref.mat

% between epochs 1 and 2
dX = XS2-XS1;
dY = YS2-YS1;
C12 = CS1+CS2;
fprintf('Epochs 1-2.\n')
fprintf('Point  dX(mm) dY(mm) \n')
for i=1:npts
    fprintf('  %2s  %6.2f %6.2f\n',Points{i},1000*dX(i),1000*dY(i));
end


%% 3rd. adjustment where reference points get identical coordinates
% (10-24a) and (10-26a)
% split of difference vector
% between epochs 1 and 2
Dr = [1000*dX(rpi);1000*dY(rpi)]; 
io = setdiff(1:npts,rpi); 
lio = length(io);
Do = [1000*dX(io);1000*dY(io)];
QS12 = QS1+QS2;
PS12 = pinv(QS12,ptol);
% split of weight matrix
irr = [rpi,rpi+npts];
ioo = [io,io+lio];
Prr = PS12(irr,irr);
Por = PS12(ioo,irr);
Poo = PS12(ioo,ioo);
Doo = pinv(Poo,ptol)*Por*Dr+Do;
QDoo = pinv(Poo,ptol);

% for object points
fprintf('Differences of object points for the 3rd adjustment\n')
fprintf('Point  dX(mm) dY(mm) sdX(mm) sdY(mm) \n')
for i=1:lio
    k = io(i);
    sdx = sqrt(C12(i,i));
    sdy = sqrt(C12(i+lio,i+lio));
    fprintf('  %2s  %6.2f %6.2f  %6.2f  %6.2f\n',Points{k},Doo(i),Doo(i+lio),sdx,sdy);
end


