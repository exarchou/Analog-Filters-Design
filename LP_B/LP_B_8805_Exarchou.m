% Dimitrios-Marios Exarchou 8805
% Designing a Low-Pass Butterwoth Filter

close all;
clear all;


%% AEM
a1 = 8;
a2 = 8;
a3 = 0;
a4 = 5;


%% Specifications
m = 1;
fp = 1.2*(3+m)*1000;
fs = 2.2*fp;

amin = 22 + (max(1,a3)-5)*(0.5);
amax = 0.65 +(max(1,a4)-3)/10;

wp = 2*pi*fp;
ws = 2*pi*fs;


%% Filter's Parameters
n = log((10^(amin/10)-1)/((10^(amax/10)-1)))/(2*log(ws/wp));
n = ceil (n);
wo = wp /((10^(amax/10)-1)^(1/(2*n)));


%% First Unit
Q1 = 0.54;
r11 = 1;
r12 = 1;
k1 = 1 + r12/r11;
C11 = 1;
R11 = 1;
R12 = Q1;
C12 = 1/Q1;
% Scaling
kf1 = wo;
C11_new = 0.1*10^(-6);
km1 = C11/(kf1*C11_new);
C11 = C11_new;
C12 = C12/(km1 * kf1);
R11 = R11*km1;
R12 = R12*km1;
r11 = r11*km1;
r12 = r12*km1;


%% Second Unit
Q2 = 1.31;
r21 = 1;
r22 = 1;
k2 = 1 + r22/r21;
C21 = 1;
R21 = 1;
R22 = Q2;
C22 = 1/Q2;
% Scaling
kf2 = wo;
C21_new = 0.1*10^(-6);
km2 = C21/(kf2*C21_new);
C21 = C21_new;
C22 = C22/(km2 * kf2);
R21 = R21*km2;
R22 = R22*km2;
r21 = r21*km2;
r22 = r22*km2;


%% Gain Adjustment
totalGain = k1*k2;
again = 10^0/totalGain;
r1 = 100;
r2 = r1*again;


%% Transfer Functions
T1 = tf([k1*wo^2], [1, wo/Q1, wo^2])
T2 = tf([k2*wo^2], [1, wo/Q2, wo^2])
TLP = series (T1 , T2);
TLP = again * TLP


%% Plotting
plot_transfer_function(T1 ,[fp wo/(2*pi) fs]);
plot_transfer_function(T2 ,[fp wo/(2*pi) fs]);
plot_transfer_function(TLP ,[fp wo/(2*pi) fs]);
plot_transfer_function(inv(TLP) ,[fp wo/(2*pi) fs]);


%% Plotting Input-Output
Fs = 100000;
t = 0:1/Fs:0.004;

T = 1/2000;
pulsewidth = 0.4*T;
pulseperiods = [0:20]*T;

u = pulstran(t,pulseperiods,@rectpuls,pulsewidth);
Vout = lsim(TLP, u, t);

figure;
plot(t,u)
axis([0 0.004 -0.5 1.5])
hold on
plot(t,Vout)
xlabel('Time(s)');
ylabel('Voltage(V)');
hold off



%% Fourier Analysis
N = length(u);
Fbins = ((0: 1/N: 1-1/N)*Fs);

U = abs(fft(u));
F = abs(fft(Vout));

helperFFT(Fbins,U)
helperFFT(Fbins,F)
