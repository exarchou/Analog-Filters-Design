% Dimitrios-Marios Exarchou 8805
% Designing a Band-Elimination Chebyshev Filter

close all;
clear all;
clc;

%% AEM
a1 = 8;
a2 = 8;
a3 = 0;
a4 = 5;


%% Specifications
f0 = 2500;
f1 = 1700 + 50*a3;
f2 = f0^2/f1;
D = (1/2.1)*(f0^2-f1^2)/f1;
f3 = (-D + sqrt(D^2 + 4*f0^2))/2;
f4 = f0^2/f3;

amin = 27 + a3*5/9;
amax = 0.4 + a4/36;


w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

Wp = 1;
Ws = (w2-w1)/(w4-w3);
bw = w2-w1;


%% Filter's Parameters
n = acosh(sqrt((10^(amin/10) - 1)/(10^(amax/10) - 1 ))) / acosh(Ws);
n = ceil(n);

e = sqrt(10^(amax/10) - 1);
a = (1/n)*(asinh(1/e));

whp = cosh(acosh(1/e)/n);


%% Butterwoth Poles
p1 = -sinh(a)*cosd(22.5)+ (1i)*cosh(a)*sind(22.5);
p2 = -sinh(a)*cosd(22.5)- (1i)*cosh(a)*sind(22.5);
p3 = -sinh(a)*cosd(67.5)+ (1i)*cosh(a)*sind(67.5);
p4 = -sinh(a)*cosd(67.5)- (1i)*cosh(a)*sind(67.5);

Wo12 = abs(p1);
Wo34 = abs(p3);
Q12 = Wo12 /(-2*real(p1));
Q34 = Wo34 /(-2*real(p3));


%% Inversing
InvW1 = 1/Wo12;
InvW2 = 1/Wo34;

P1 = real(p1)/Wo12^2 + (1i)*imag(p1)/Wo12^2;
P2 = real(p1)/Wo12^2 - (1i)*imag(p1)/Wo12^2;
P3 = real(p3)/Wo34^2 + (1i)*imag(p3)/Wo34^2;
P4 = real(p3)/Wo34^2 - (1i)*imag(p3)/Wo34^2;


%% Geffe Algorithm
qc = w0/bw;

% Pole 1 transformation
S12 = -real(P1);
W12 = imag(P1);

C1 = S12^2 + W12^2;
D1 = 2*S12/qc;
E1 = 4 + C1/qc^2;
G1 = sqrt(E1^2 - 4*D1^2);
Q_12 = sqrt((E1 + G1)/2)/D1;
k1 = (S12*Q_12)/qc;
W1 = k1 + sqrt(k1^2 - 1);
wo1 = w0/W1;
wo2 = W1*w0;

% Pole 3 transformation
S34 = -real(P3);
W34 = imag(P3);

C2 = S34^2 + W34^2;
D2 = 2*S34/qc;
E2 = 4 + C2/qc^2;
G2 = sqrt(E2^2 - 4*D2^2);
Q_34 = sqrt((E2 + G2)/2)/D2;
k2 = (S34*Q_34)/qc;
W2 = k2 + sqrt(k2^2 - 1);
wo3 = w0/W2;
wo4 = W2*w0;


%% Zeros
wz1 = w0;
wz2 = w0;
wz3 = w0;
wz4 = w0;


%% First Unit LPN Boctor
wzo1 = wz1/wo1;
kmin1 = wo1^2/wz1^2;
k11 = 0.8;

R11 = 2/(k11*wzo1^2 - 1);
R12 = 1/(1 - k11);
R13 = (k11/Q_12^2 + k11*wzo1^2 - 1)/2;
R14 = 1/k11;
R15 = 1;
R16 = 1;
C11 = k11/(2*Q_12);
C12 = 2*Q_12;
k1 = 1/((1/2)*(k11/Q_12^2 + k11*wzo1^2 + 1));
 
% Scaling
kf1 = wo1;
km1 = C11/(kf1 * 0.01*10^(-6));
R11 = R11*km1;
R12 = R12*km1;
R13 = R13*km1;
R14 = R14*km1;
R15 = R15*km1;
R16 = R16*km1;
C11 = 0.01*10^(-6);
C12 = C12/(km1*kf1);


%% Second Unit HPN Boctor
wzo2 = wz2/wo2;
BoctorHighPass( wz2, wo2, Q_12);
k2 = circuit.H;


%% Third Unit LPN Boctor
wzo3 = wz3/wo3;
kmin3 = wo3^2/wz3^2;
k31 = 0.7;

R31 = 2/(k31*wzo3^2 - 1);
R32 = 1/(1 - k31);
R33 = (k31/Q_34^2 + k31*wzo3^2 - 1)/2;
R34 = 1/k31;
R35 = 1;
R36 = 1;
C31 = k31/(2*Q_34);
C32 = 2*Q_34;
k3 = 1/((1/2)*(k31/Q_34^2 + k31*wzo3^2 + 1));

% Scaling
kf3 = wo3;
km3 = C31/(kf3 * 0.01*10^(-6));
R31 = R31*km3;
R32 = R32*km3;
R33 = R33*km3;
R34 = R34*km3;
R35 = R35*km3;
R36 = R36*km3;
C31 = 0.01*10^(-6);
C32 = C32/(km3*kf3);


%% Fourth Unit High Pass Notch 7.21
%BoctorHighPass( wz4, wo4, Q_34)  
wzo4 = wz4/wo4;
k41 = 1/wzo4^2 - 1;
k42 = (2+k41)*Q_34^2/((2+k41)*Q_34^2 + 1);

R41 = 1;
R42 = (2+k41)^2 * Q_34^2;
R43 = 1;
R44 = (2+k41) * Q_34^2;
C42 = 1/((2+k41)*Q_34);
C41 = k41*C42;
k4 = k42/wzo4^2;

% Scaling
kf4 = wo4;
km4 = C42/(kf4* 0.01*10^(-6));
R41 = R41*km4;
R42 = R42*km4;
R43 = R43*km4;
R44 = R44*km4;
C41 = C41/(km4*kf4);
C42 = 0.01*10^(-6);



%% Transfer Functions
T1 = tf( [k1 0 (k1*wz1^2)], [ 1 (wo1/Q_12) wo1^2 ] )
T2 = tf( [k2 0 (k2*wz2^2)], [ 1 (wo2/Q_12) wo2^2 ] )
T3 = tf( [k3 0 (k3*wz3^2)], [ 1 (wo3/Q_34) wo3^2 ] )
T4 = tf( [k4 0 (k4*wz4^2)], [ 1 (wo4/Q_34) wo4^2 ] )


%Overall TF
TBE = series(T1, T2);
TBE = series(TBE, T3);
TBE = series(TBE, T4);


%% Gain Adjustment at 5 dB
totalGain = abs(evalfr(TBE, 0*1i));
again = (10^(0.25))/totalGain;
TBE = again * TBE
Z1 = 10^4;
Z2 = Z1*again;


%% Plotting
plot_transfer_function( T1, [f1 f3 f0 f4 f2] );
plot_transfer_function( T2, [f1 f3 f0 f4 f2] );
plot_transfer_function( T3, [f1 f3 f0 f4 f2] );
plot_transfer_function( T4, [f1 f3 f0 f4 f2] );
plot_transfer_function( TBE, [f1 f3 f0 f4 f2] );
plot_transfer_function( inv(TBE), [f1 f3 f0 f4 f2] );


%% Plotting Input-Output
fa = f0 - (f0-f3)/2;
fb = f0 + (f0+f3)/2;
fc = 0.5*f1;
fd = 2.4*f2;
fe = 3.5*f2;

Fs = 100000;
t = 0:1/Fs:0.003;

u = 0.8*cos(2*pi*fa*t)+cos(2*pi*fb*t)+cos(2*pi*fc*t)+0.8*cos(2*pi*fd*t)+0.4*cos(2*pi*fe*t);
Vout = lsim(TBE, u, t);

figure;
plot(t,u)
hold on
plot(t,Vout)
xlabel('Time(s)');
ylabel('Voltage(V)');
hold off


%% Fourier Analysis
N = length(u);
Fbins = (0: 1/N: 1-1/N)*Fs;

U = abs(fft(u));
F = abs(fft(Vout));

helperFFT(Fbins,U)
axis([0 5*10^4 0 200])
helperFFT(Fbins,F)
axis([0 5*10^4 0 200])
