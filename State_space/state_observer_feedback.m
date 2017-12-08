%% 1.2.2.3. Example Implementation [Ostertag]
% State observer + state feedback control
% Date: 24.09.2016

clear all
clc

%% Continuous plant model
Ac = [0, 1, 0; 0, -1, 1; 0, 0, -2];
Bc = [0; 0; 1];
C = [1, 0, 0];
D = 0;

%% Controlability matrix
CO = ctrb(Ac, Bc);
display(rank(CO), 'rank[Co(Ac,Bc)]');

%% Observability matrix
OB = obsv(Ac, C);
display(rank(OB), 'rank[Ob(Ac,C)]');

%% Plant discretization
% Tustin's method

Ts = 0.1; % sample time
I = eye(3);

% first order Pade approximation
Ad = (I-Ac*Ts/2)\(I+Ac*Ts/2);   
Bd = (I-Ac*Ts/2)\Bc*Ts;
% C, D do not change after discretization

%% Luenberger observer
L = acker(Ad', C', [0.45, 0.45, 0.5])';

%% State feedback for the controllable system
display(eig(Ad)', 'Open-loop system eigenvalues');

% desired closed-loop poles
p = [0.5, 0.5+0.2*1i, 0.5-0.2*1i];
% use Ackermann's formula
K = acker(Ad, Bd, p); % for transient response; regulation behaviour

A0 = Ad-Bd*K;
display(eig(A0)', 'Closed-loop system eigenvalues');

% for steady-state response; servo behaviour
% H(s)|s=0 = 1
% H(z)|z=1 = 1
% M = inv(C*inv(I-Ad+Bd*K)*Bd); % static gain
M = inv(C/(I-Ad+Bd*K)*Bd);

%% Testing
x = [5; -1; 1];
xhat = [0; 0; 0];

XX = x;
XXhat = xhat;

N = 60; % number of samples to consider
YR = [1.0*ones(1,N/2+1), 3.0*ones(1,N/2)]; % input signal; reference

for k=0:N-1  % "current" estimator
    % (k+1) - present sample
    % ( k ) - previous sample
    
    % y computation placement has implications with respect to non-zero
    % computation time
    yr = YR(k+1);
    y = C*x+D*yr;       % y(k+1) = C*x(k)+D*yr(k+1)
    yhat = C*xhat+D*yr; % yhat(k+1) = C*xhat(k)+D*yr(k)
    
    % u - separated to avoid instability from saturation
    u = -K*xhat+M*yr;   % u(k) = -K*xhat(k)+M*yr(k+1)
    x = Ad*x+Bd*u;      % x(k+1) = A*x(k)+B*u(k)
    % xhat(k+1) = A*xhat(k)+B*u(k)+L*(y(k+1)-yhat(k+1)
    xhat = Ad*xhat+Bd*u+L*(y-yhat); 
    
    XX = [XX, x];
    XXhat = [XXhat, xhat];
end

plot((0:N), [YR; XX(1,:); XXhat(1,:)]);
title('Model-in-the-Loop Simulation');
xlabel('n [Samples]'); ylabel('Output');
legend('Set-Point', 'Real', 'Observer', 'location', 'southeast')