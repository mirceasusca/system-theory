%% Simple modal control of a continuous plant
% No delay
% Completely controllable and observable system
% Eric Ostertag 
%
% n - number of states
% p - number of inputs
% q - number of outputs
%
% Date: 15.10.2016

clear all
close all
clc

%% Non-canonical form
Ac=[-2.8428   -3.6378    4.1648;
    0.6879   -6.5990    2.9545;
    0.0459    3.8149   -6.1913];
%
Bc=[-0.4622, 0.234;
    0.1772,  -0.1189;
    1.1079, 0.82];  
% Bc=[-0.4622;
%     0.1772;
%     1.1079]; 
C=[1, 1, 0.5];
%
D=[0, 0];

%% Plant discretization
% Tustin's method
Ts = 0.1; % sample time
I = eye(3);
% 
% first order Pade approximation
Ad = (I-Ac*Ts/2)\(I+Ac*Ts/2);   
Bd = (I-Ac*Ts/2)\Bc*Ts;
% C, D do not change after discretization

%% Diagonal cannonical form; starting point at 1.4.2, pg. 23
[T,Diag]=eigs(Ad,size(Ad,1));   % sorted ascendingly
% shift the poles closer to the imaginary axis (big settling-time)
Diag=Diag(:,1:length(Diag));
Diag=Diag(1:length(Diag),:);    % sorted descendingly
%
AA=T\Ad*T;
BB=T\Bd;
CC=C*T;
DD=D;
%
[~,n] = size(BB);
p = n;                          % number of inputs
clear n
% 
[m,~] = size(CC);               % number of outputs
q = m;
clear m
%  
[m,n] = size(AA);               % number of states
if m ~= n
   error('The state matrix must be square.'); 
end
clear m

%% Model-in-the-Loop (MiL) simulation: open-loop system
% % Variable declarations and algorithm
% N = 26;                         % number of simulation points
% % 
% x = [0; 0; 0];                  % initial state
% % 
% t = (0:N-1)*Ts;
% u = [0*ones(1,length(t));       % inputs
%      ones(1,length(t))];
% 
% XX = zeros(n,N);
% XX(:,1) = x;
% YY = zeros(q,N);
% 
% for i=1:N,
%    % x[k+1]=A*x[k]+B*u[k]
%    % y[k]  =C*x[k]+D*u[k]
%     
%    y = CC*x+DD*u(:,i);
%    x = AA*x+BB*u(:,i);
% 
%    XX(:,i) = x;
%    YY(:,i) = y;
% end

%% Plots
% plot_states(t,XX);
% plot_outputs(t,YY);

%% Modal controller design
% Split in controlled part and uninfluenced part
% p -> number of controlled modes (same as rank of the input matrix BB)
% n-p -> remaining uncontrolled modes
Bhatp = BB(1:p,:);              % size pxp
Diagp = Diag(1:p,1:p);          % size pxp
DiagpRef = diag([0.42, 0.46]);
%
Tinv = inv(T);
Tinvp = Tinv(1:p,:);        % first p rows of the inverse transf. matrix
%
L = (Tinvp*Bd)\(Diagp-DiagpRef)*Tinvp;  % Simple modal controller (pxn)

%% Modal control, (MiL) simulation: closed-loop system
% Variable declarations and algorithm
N = 51;                         % number of simulation points
% 
% System initialization
x = [2.1; -0.20; 0.7];                  % initial state
%
XX = zeros(n,N);
YY = zeros(q,N);
XX(:,1) = x;
%
t = (0:N-1)*Ts;
YR = [ones(1,length(t)-27), 0.83*ones(1,27);        % reference inputs
      ones(1,length(t))];
%
for k=2:N
    yref = YR(:,k);
    y=C*x+D*yref;

    % Regulation behaviour
    u=-L*x;  
    
    % Servo behaviour
    if p<q
        % TODO cannot servo all outputs
        % drop columns from C
    elseif p==q
        M=inv(CC/(I+Bd*L-Ad)*Bd);
        u=u+M*yref;
    elseif p>q
        % drop p-q control inputs for the servo behaviour
        Bdsrv=Bd(:,1:p-q);
        M=inv(C/(I+Bd*L-Ad)*Bdsrv);
        u(1:p-q)=u(1:p-q)+M*yref(1:p-q);
    end
    
    % system eigenvalues are from matrix (AA-BB*L)
    x=Ad*x+Bd*u;
    
    XX(:,k) = x;
    YY(:,k) = y;
end

%% Plots
plot_states(t,XX);
plot_outputs(t,YY);