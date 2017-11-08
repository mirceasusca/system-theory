%% Sallen-Key Notch Filter
close all
clear all
clc

%%
R1 = 3.3e+03;
R2 = 3.2e+03;
R3 = 1.65e+03;
R4 = 1e+02;
R5 = 0; % =0; eliminates the negative voltage divider feedback; Q = 1/2
C1 = 100e-09;
C2 = 100e-09;
C3 = 200e-09; 

%% Fixed Q=1/2; Voltage follower configuration
% A = [-1/R3/C1, 0, -1/R2/C1;
%      0, 0, -1/R2/C2;
%      1/R1/C3, 1/R1/C3, -(1/R1+1/R2)/C3];
% %
% B = [1/R3/C1;
%     0;
%     0];
% %
% C = [-1,-1,0];
% %
% D = 1;   

%% Variable Q; Voltage divider on inverting input
A = [(R5/R4/R2-1/R3)/C1, R5/R4/R2/C1, -1/R2/C1;
     R5/R4/R2/C2, R5/R4/R2/C2, -1/R2/C2;
     (1/R1+R5/R4*(1/R1+1/R2))/C3, (1/R1+R5/R4*(1/R1+1/R2))/C3, -(1/R1+1/R2)/C3];
%
B = [(1/R3-R5/R4/R2)/C1;
    -R5/R4/R2/C2;
    -R5/R4*(1/R1+1/R2)/C3];
%
C = (1+R5/R4)*[-1,-1,0];
%
D = (1+R5/R4);   

%% Frequency response of the third-order filter
f = logspace(2,3,16000);
sys = ss(A,B,C,D);
[mag,ph] = bode(sys,f*2*pi);

%% Ideal second-order filter comparison
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

subplot(212);
semilogx(f,ph(:),'k'); grid; hold on;
semilogx(f,(ph2(:)),'r'); 
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');