t = ElectricalSignals_cl.time;

irv_cl = ElectricalSignals_cl.signals(1).values;
isv_cl = ElectricalSignals_cl.signals(2).values;
phi_s_v_cl = ElectricalSignals_cl.signals(3).values;
vs_cl = ElectricalSignals_cl.signals(4).values;
wm_cl=MechanicalSignals_cl.signals(1).values;
Te_cl=MechanicalSignals_cl.signals(2).values;

irv_cli = ElectricalSignals_cli.signals(1).values;
isv_cli = ElectricalSignals_cli.signals(2).values;
phi_s_v_cli = ElectricalSignals_cli.signals(3).values;
vs_cli = ElectricalSignals_cli.signals(4).values;
wm_cli=MechanicalSignals_cli.signals(1).values;
Te_cli=MechanicalSignals_cli.signals(2).values;

irv_ol = ElectricalSignals_ol.signals(1).values;
isv_ol = ElectricalSignals_ol.signals(2).values;
phi_s_v_ol = ElectricalSignals_ol.signals(3).values;
vs_ol = ElectricalSignals_ol.signals(4).values;
wm_ol=MechanicalSignals_ol.signals(1).values;
Te_ol=MechanicalSignals_ol.signals(2).values;
T_L = TorqueData(:,2);
%%
figure;
ax1=subplot(426);plot(t,irv_cl);
title('SVM PWM input rotor currents')

ax3=subplot(424);plot(t,isv_cl);
title('SVM PWM input stator currents')

ax5=subplot(428);plot(t,phi_s_v_cl);
title('SVM PWM input stator flux signals')
xlabel('Time [s]');

ax7=subplot(422);plot(t,vs_cl)
title('(b) SVM PWM input stator voltages')
%
ax2=subplot(425);plot(t,irv_ol);
title('Sine input rotor currents')
ylabel('i_r abc [A]');

ax4=subplot(423);plot(t,isv_ol);
title('Sine input stator currents')
ylabel('i_s abc [A]');

ax6=subplot(427);plot(t,phi_s_v_ol);
title('Sine input stator flux signals')
ylabel('\Phi abc [Wb]');
xlabel('Time [s]');

ax8=subplot(421);plot(t,vs_ol)
title('(a) Sine input stator voltages')
ylabel('V_s abc [V]');
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'x')

%%
figure;
ax1=subplot(211);plot(t,wm_cl(:,1),'--');hold
plot(t,wm_ol(:,2),'linewidth',1.2)
plot(t,wm_cl(:,2),':','linewidth',1.5);
plot(t, wm_cli(:,2),'k','linewidth',1.0);
legend('\omega_m^*', '\omega_m^{OL}','\omega_m^{F}','\omega_m^{C}');
ylabel('\omega_m [rpm]')
xlabel('Time [s]'); grid minor
title('Mechanical subsystem');

ax2=subplot(212);plot(t,T_L,'--',t,Te_ol,t,Te_cl, t, Te_cli);
legend('T_L','T_e^{OL}','T_e^{F}','T_e^{C}');
ylabel('Electromechanical torque [Nm]');
xlabel('Time [s]');
linkaxes([ax1,ax2],'x')

%%
figure;
ax1=subplot(311);plot(t,vs_cl(:,1));
ax2=subplot(312);plot(t,vs_cl(:,3));
ax3=subplot(313);plot(t,vs_cl(:,2));
linkaxes([ax1,ax2,ax3],'x')

%%
figure;
t=SPT_data.time;
speed=SPT_data.signals(1).values;
power=SPT_data.signals(2).values;
torque=SPT_data.signals(3).values;
% (15000:end)

ax1=subplot(211);plot(t(15000:end),speed(15000:end));
ax2=subplot(212);
yyaxis right;
plot(t(15000:end),power(15000:end));
yyaxis left;
plot(t(15000:end),torque(15000:end)); grid minor;
linkaxes([ax1,ax2],'x')

%% SVPWM
t=SVPWM_wave.time;
mod_wave=SVPWM_wave.signals(1).values;
transistor_comm=SVPWM_wave.signals(2).values;
l2l_voltage=SVPWM_wave.signals(3).values;

figure;
ax1=subplot(311); plot(t,mod_wave(:,2),'-','linewidth',0.25); ylim([-1.2,1.2]); hold
plot(t,mod_wave(:,1))
title('Space Vector PWM Implementation'); ylabel('Modulation wave');
ax2=subplot(312); plot(t,transistor_comm); ylim([-0.2,1.2])
ylabel('Commutation wave');
ax3=subplot(313); 
plot(t, 220*sqrt(2)*l2l_voltage); hold on
plot(t, 220*sqrt(2)*sin(2*pi*40*t+pi/1.1965),'linewidth',1.0)
ylabel('Line-to-line voltage');
% plot(t, 220*sqrt(2)*sin(2*pi*40*t+pi/1.23))
% plot(t, 220*sqrt(2)*sin(2*pi*40*t+pi/1.25))
xlabel('Time [s]');
% plot(t, 220*sqrt(2)*sin(2*pi*40*t+pi/1.27))
% legend('1','2','3','4','5')
ylim([-250*sqrt(2), 250*sqrt(2)])
linkaxes([ax1,ax2,ax3],'x')

%% 
figure;
t = 0:0.1:100;
x = 3.5*t;

for i=1:length(t)
%    if x(i) <=10
%        x(i) = 10;
%    end
   if t(i) < 10
       x(i)=10*3.5;
   end
   if x(i) >=220*sqrt(2)
       x(i)=220*sqrt(2);
   end
end

plot(0:0.1:20,3.5*[0:0.1:20],'--','linewidth',1.5);hold on
plot(t,x,'linewidth',1.2)
title('Constant V/Hz characteristic')
xlabel('Frequency [Hz]');
ylabel('Voltage [V]'); grid minor
% txt1 = '\leftarrow sin(\pi) = 0';
text(1,4.6*12,'f_{min}= 10')
text(84,95*3.5,'V_{max}= 311')