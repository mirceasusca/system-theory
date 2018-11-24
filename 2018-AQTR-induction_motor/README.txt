  The model has been developed using MATLAB 2017a. 
  If you download the model using different (older) versions of MATLAB, 
there may be different settings for some of the blocks.

  In case this happens, to reproduce the results from the paper, we suggest the following settings:
  1) The Solver has been set as Fixed-Step with a Fundamental Step Time of Ts = 1e-5 seconds
for this example. The PowerGUI block also is Discrete with a sampling-time of 1e-5 s. The integration
method has been set to Fixed-step with the Solver selected as auto (Automatic solver selection), but 
ode3 (Bogacki-Shampine) works well.
  2) The Repeating Sequence block for the SVPWM has the Time Values [0 0.25 0.5 0.75 1]/(10e3)
and Output Values [1 0 -1 0 1], as seen in "repeating_sequence_settings.png".
  3) The Vdc block from the MosFET bridge has the voltage set at 1 V, as seen in "vdc_settings_1v.png".
  4) Final time set to at least 0.4 seconds.