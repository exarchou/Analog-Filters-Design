% Dimitrios-Marios Exarchou 8805
% Designing a High-Pass Chebyshev Filter

close all;
clear all;


%% AEM
a1 = 8;
a2 = 8;
a3 = 0;
a4 = 5; 


%% Specifications
m = 2;
fp = (3+m)*10^3;
fs = fp/1.8;

amin = 25 + a3*4/9;
amax = 0.5 + a4*0.25/9;

wp = 2*pi*fp;
ws = 2*pi*fs;


%% Transformation
Wp = 1;
Ws = wp/ws;


%% Filter's Parameters
e = sqrt(10^(amax/10)-1);
n = (acosh(((10^(amin/10)-1)/(10^(amax/10)-1))^(1/2)))/acosh(Ws);
n = ceil(n);
a = (asinh(1/e))/n;
Whp = cosh((acosh((10^(amax/10)-1)^(-1/2)))/n);


%% Butterwoth Poles
p1 = -sinh(a)*cosd(22.5)+ (1i)*cosh(a)*sind(22.5);
p2 = -sinh(a)*cosd(22.5)- (1i)*cosh(a)*sind(22.5);
p3 = -sinh(a)*cosd(67.5)+ (1i)*cosh(a)*sind(67.5);
p4 = -sinh(a)*cosd(67.5)- (1i)*cosh(a)*sind(67.5);

Wo12 = sqrt((real(p1))^2+(imag(p1))^2);
Q12 = Wo12 /(2 * abs(real(p1)));

Wo34 = sqrt((real(p3))^2+(imag(p3))^2);
Q34 = Wo34/(2*abs(real(p3)));


%% Inversing
w_hp = wp/Whp;
w12 = wp/Wo12;
w34 = wp/Wo34;


%% First Unit
C11 = 1;
C12 = C11;
R11 = 1;
R12 = R11;
k1 = 3-1/Q12;
r11 = 1;
r12 = 2-1/Q12;
% Scaling
kf1 = w12;
C11_new = 1*10^(-6);
km1 = C11/(C11_new*kf1);

R11 = R11 * km1;
R12 = R12 * km1;
C12 = C12/(kf1*km1);
C11 = C11_new;
r11 = r11 * km1;
r12 = r12 * km1;


%% Second Unit
C21 = 1;
C22 = C21;
R21 = 1;
R22 = R21;
k2 = 3-1/Q34;
r21 = 1;
r22 = 2-1/Q34;
% Scaling
kf2 = w34;
C21_new = 1*10^(-6);
km2 = C21/(C21_new*kf2);

R21 = R21 * km2;
R22 = R22 * km2;
C22 = C22/(kf2*km2);
C21 = C21_new;
r21 = r21 * km2;
r22 = r22 * km2;


%% Transfer Functions
T1 = tf ([ k1 0 0], [1 w12/Q12 w12^2])
T2 = tf ([ k2 0 0], [1 w34/Q34 w34^2])
THP = T1*T2;


%% Gain Adjustment
totalGain = k1*k2;
again = (10^(0))/totalGain; 
THP = again * THP
r1 = 100;
r2 = r1*again;


%% Plotting
plot_transfer_function(T1, [fs w_hp/(2*pi) fp]);
plot_transfer_function(T2, [fs w_hp/(2*pi) fp]);
plot_transfer_function(THP, [fs w_hp/(2*pi) fp]);
plot_transfer_function(inv(THP), [fs w_hp/(2*pi) fp]);


%% Plotting Input-Output
f1 = 0.2*ws/(2*pi);
f2 = 0.7*ws/(2*pi);
f3 = 1.6*wp/(2*pi);
f4 = 2.4*wp/(2*pi);
f5 = 5*wp/(2*pi);

Fs = 100000;
t = 0:1/Fs:0.002;

u = cos(2*pi*f1*t)+0.6*cos(2*pi*f2*t)+1.5*cos(2*pi*f3*t)+0.7*cos(2*pi*f4*t)+0.4*cos(2*pi*f5*t);
Vout = lsim(THP, u, t);

figure;
plot(t,u)
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
axis([0 5*10^4 0 160])
helperFFT(Fbins,F)
axis([0 5*10^4 0 160])
