%% Sallen-Key Lowpass filter with fixed Q
R1 = 3.3e+03;
R2 = 3.3e+03;
C1 = 100e-09;
C2 = 100e-09;

%% Fixed Q; voltage follower configuration
A = [0, -1/R2/C1;
     1/R1/C2, -(1/R1+1/R2)/C2];
%
B = [0;
     -1/R1/C2];
% 
C = [1, 0];
%
D = 0;

%%
f = logspace(1,4,1000);
sys = ss(A,B,C,D);
[mag,ph] = bode(sys,2*pi*f);

figure;
%
subplot(211); semilogx(f,20*log10(mag(:))); grid; hold on;
semilogx(f,(max(20*log10(mag(:)))-3)*ones(1,length(f)),'-r');
title('High-pass Sallen-Key filter with variable gain');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
%
subplot(212); semilogx(f,ph(:)); grid;
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');