%% Dual-operational amplifier Twin-T Notch Filter
close all; clear all; clc

%%
R1 = 3.3e+03;
R2 = 3.2e+03;
R3 = 1.65e+03;
C1 = 100e-09;
C2 = 100e-09;
C3 = 200e-09;
%
% R5 = R4 => Q = 1/2
R4 = 100;
R5 = 100; 

%% Variable Q
A = [-R4/(R4+R5)*(1/R2+1/R3)/C1, (-1/R2+R5/(R4+R5)*(1/R2+1/R3))/C1, 1/R2/C1;
    -R4/(R4+R5)/R2/C2, -R4/(R4+R5)/R2/C2, 1/R2/C2;
    (1/R2-R5/(R4+R5)*(1/R1+1/R2))/C3, (1/R2-R5/(R4+R5)*(1/R1+1/R2))/C3, -(1/R1+1/R2)/C3];
%
B = [R4/(R4+R5)*(1/R2+1/R3)/C1;
    R4/(R4+R5)/R2/C2;
    -R4/(R4+R5)*(1/R1+1/R2)/C3];
% 
C = [-1, -1, 0];
%
D = 1;

%% Frequency response
f = logspace(1,4,16000);
sys = ss(A,B,C,D);
[mag,ph] = bode(sys,f*2*pi);

%% Comparison with ideal second-order notch filter
R = R1;
C = C1;
H2 = tf([1,0,(1/R/C)^2],[1,2/R/C,(1/R/C)^2]);
[m2,ph2] = bode(H2,f*2*pi);

figure;
subplot(211);
semilogx(f,20*log10(m2(:)),'r');  grid; hold on;
semilogx(f,20*log10(mag(:)),'k');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('Bode diagram');
legend('H_2(s)','H_3(s)','Location','best');
%
subplot(212);
semilogx(f,(ph2(:)),'r'); grid; hold on;
semilogx(f,ph(:),'k');
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');