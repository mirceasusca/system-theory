%% Motor parameters (Bodine Electric Company model 295)
clear all;
%
P=4;
Rs=14.6;
Lls=0.0222;
Rr=12.76;
Llr=0.0518;
Lm=0.2963;
J=0.001;
B=0.000124;
%
Lmlstar=1/(1/Lls+1/Lm+1/Llr);
