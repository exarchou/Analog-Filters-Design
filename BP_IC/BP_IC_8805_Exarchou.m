% Dimitrios-Marios Exarchou 8805
% Designing a Band-Pass Inverse Chebyshev Filter

close all;
clear all;


%% AEM
a1 = 8;
a2 = 8;
a3 = 0;
a4 = 5;


%% Specifications
f0 = 1650;
f1 = 1400 + 25*a4;
f2 = f0^2/f1;
D = 2.5*(f0^2-f1^2)/f1;
f3 = (-D + sqrt(D^2 + 4*f0^2))/2;
f4 = f0^2/f3;

amin = 34 + a3*5/9;
amax = 0.4 + a4/36;

w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

Wp = 1;
Ws = (w4-w3)/(w2-w1);
bw = w2-w1;


%% Filter's Parameters
n = acosh(sqrt((10^(amin/10) - 1)/(10^(amax/10) - 1 ))) / acosh(Ws);
n = ceil(n);

e = 1/sqrt(10^(amin/10) - 1);
a = (1/n)*(asinh(1/e));

whp = 1/cosh(acosh(1/e)/n);


%% Poles
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

InvW1 = InvW1*Ws;
InvW2 = InvW2*Ws;

S12 = -InvW1/(2*Q12);
S34 = -InvW2/(2*Q34);
W12 = sqrt(InvW1^2 - S12^2);
W34 = sqrt(InvW2^2 - S34^2);


%% Zeros 
Z1 = sec(pi/8);
Z2 = sec(3*pi/8);

Z1 = Z1*Ws;
Z2 = Z2*Ws;


%% Geffe Algorithm
qc = w0/bw;

% Pole 1 transformation
C1 = S12^2 + W12^2;
D1 = -2*S12/qc;
E1 = 4 + C1/qc^2;
G1 = sqrt(E1^2 - 4*D1^2);
Q_12 = sqrt((E1 + G1)/2)/D1;
k1 = (-S12*Q_12)/qc;
W1 = k1 + sqrt(k1^2 - 1);
wo1 = w0/W1;
wo2 = W1*w0;

% Pole 3 transformation
C2 = S34^2 + W34^2;
D2 = -2*S34/qc;
E2 = 4 + C2/qc^2;
G2 = sqrt(E2^2 - 4*D2^2);
Q_34 = sqrt((E2 + G2)/2)/D2;
k2 = (-S34*Q_34)/qc;
W2 = k2 + sqrt(k2^2 - 1);
wo3 = w0/W2;
wo4 = W2*w0;


% Zero 1 transformation
Kzero1 = 2 + (Z1^2)/(qc^2);
x1 = (Kzero1 + sqrt(Kzero1^2 - 4))/2;
wz1 = w0*sqrt(x1);
wz2 = w0/sqrt(x1);

% Zero 2 transformation
Kzero2 = 2 + (Z2^2)/(qc^2);
x2 = (Kzero2 + sqrt(Kzero2^2 - 4))/2;
wz3 = w0*sqrt(x2);
wz4 = w0/sqrt(x2);



%% First Unit LPN 
wzo1 = wz1/wo1;
R11 = 1;
R12 = 4*Q_12^2;
R13 = wzo1^2 / (2*Q_12^2);
R14= 1;
R15= 4*Q_12^2/(wzo1^2 - 1);
C1 = 1/(2*Q_12);
k1 = 1/(1 +  wzo1^2/(2*Q_12^2));
% Scaling
kf1 = wo1;
km1 = C1/(kf1 * 10^(-6));
R11 = R11*km1;
R12 = R12*km1;
R13 = R13*km1;
R14 = R14*km1;
R15 = R15*km1;
C1 = 10^(-6);


%% Second Unit HPN 
wzo2 = wz2/wo2;
k21 = 1/wzo2^2 - 1;
k22 = (2 + k21)*Q_12^2/((2 + k21)*Q_12^2 + 1);
R21 = 1;
R22 = (2 + k21)^2*Q_12^2;
R23 = 1;
R24 = (2 + k21)*Q_12^2;
C22 = 1/((2+k21) * Q_12);
C21 = k21*C22;
k2 = k22/wzo2^2;
% Scaling
kf2 = wo2;
km2 = C22/(kf2*10^(-6));
R21 = R21*km2;
R22 = R22*km2;
R23 = R23*km2;
R24 = R24*km2;
C21 = C21/(km2*kf2);
C22 = 10^(-6);


%% Third Unit LPN 
wzo3= wz3/wo3;
R31 = 1;
R32 = 4*Q_34^2;
R33 = wzo3^2/(2*Q_34^2);
R34 = 1;
R35 = 4*Q_34^2/(wzo3^2 - 1);
C3 = 1/(2*Q_34);
k3 = 1/( 1 + wzo3^2/(4*Q_34^2));
% Scaling
kf3 = wo3;
km3 = C3/(kf3*10^(-6));
R31 = R31*km3;
R32 = R32*km3;
R33 = R33*km3;
R34 = R34*km3;
R35 = R35*km3;
C3 = 10^(-6);


%% Fourth Unit HPN 
wzo4 = wz4/wo4;
k41 = 1/wzo4^2 - 1;
k42 = (2 + k41)*Q_34^2/((2 + k41)*Q_34^2 + 1);
R41 = 1;
R42 = (2 + k41)^2*Q_34^2;
R43 = 1;
R44 = (2 + k41)*Q_34^2;
C42 = 1/((2+k41)*Q_34);
C41 = k41*C42;
k4 = k42/wzo4^2;
% Scaling
kf4 = wo4;
km4 = C42/(kf4*10^(-6));
R41 = R41*km4;
R42 = R42*km4;
R43 = R43*km4;
R44 = R44*km4;
C41 = C41/(km4*kf4);
C42 = 10^(-6);


%% Transfer Functions
T1 = tf( [k1 0 ( k1 * wz1^2 ) ], [ 1 ( wo1 / Q_12 ) wo1^2 ] )
T2 = tf( [k2 0 ( k2 * wz2^2 ) ], [ 1 ( wo2 / Q_12 ) wo2^2 ] )
T3 = tf( [k3 0 ( k3 * wz3^2 ) ], [ 1 ( wo3 / Q_34 ) wo3^2 ] )
T4 = tf( [k4 0 ( k4 * wz4^2 ) ], [ 1 ( wo4 / Q_34 ) wo4^2 ] )


%Overall TF
TBP = series(T1, T2);
TBP = series(TBP, T3);
TBP = series(TBP, T4);


%% Gain Adjustment at 5 dB
totalGain = abs(evalfr(TBP, w0*1i));
again = (10^(0.25))/totalGain;
TBP = again * TBP
r1 = 10^4;
r2 = r1*again;


%% Plotting
plot_transfer_function( T1, [f3 f1 f0 f2 f4] )
plot_transfer_function( T2, [f3 f1 f0 f2 f4] )
plot_transfer_function( T3, [f3 f1 f0 f2 f4] )
plot_transfer_function( T4, [f3 f1 f0 f2 f4] )
plot_transfer_function( TBP, [f3 f1 f0 f2 f4] )
plot_transfer_function( inv(TBP), [f3 f1 f0 f2 f4] )



%% Plotting Input-Output
fa = f0 - (f0-f1)/3;
fb = f0 + (f0+f1)/4;
fc = 0.5*f3;
fd = 2.4*f4;
fe = 3*f4;

Fs = 100000;
t = 0:1/Fs:0.005;

u = cos(2*pi*fa*t)+0.6*cos(2*pi*fb*t)+0.7*cos(2*pi*fc*t)+0.8*cos(2*pi*fd*t)+0.6*cos(2*pi*fe*t);
Vout = lsim(TBP, u, t);

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
axis([0 3*10^4 0 300])
helperFFT(Fbins,F)
axis([0 3*10^4 0 300])