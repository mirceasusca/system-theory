%% Script which transforms a system to controllability canonical form
% dx/dt = Ax+bu => dz/dt = Ac*z+bc*u
%
% V = [qc';
%     qc'*A;
%     qc'*A^2;
%     ...
%     qc'*A^n-1]
%
% qc' is the last row of the inverse controllability matrix Qc^-1, deduced 
% from the following system of equations:
%   qc'*b=0
%   qc'*A*b=0
%   ...
%   qc'*A^(n-1)*b=0
% [also: qc'*Qc=(0 ... 0 1)]
%
% Date: 09.10.2016

clear all
close all
clc

%% Non-canonical form
A=[-0.6980    3.8213    2.1263;
    0.0526   -0.5982   -0.0151;
   -0.8929    0.6250    0.0961];

b=[-0.4622;
    0.1772;
   -0.1079];

%%
Qc = ctrb(A,b);
QcInv = inv(Qc); % last row
qc = QcInv(length(QcInv),:)';

V = [qc';
    qc'*A;
    qc'*A^2];

% Controllability canonical form
Ac = V*A/V;
bc = V*b;