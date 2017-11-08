function [RH, is_stable] = rh(P, varargin)
% RH - Calculates the Routh-Hurwitz matrix for the characteristic polynomial P.
% Calculates the number of poles in the right-half plane and evaluates the
% stability of the closed-loop system.
%
% Syntax:
%     rh(P,eps);
%     rh(P);
%     RH = rh(P,eps)
%     RH = rh(P);
%
% Inputs:
%    P   - coefficient vector in the decreasing order of the polynomial powers
%    eps - value used when there's a zero on the first column (considered as a perturbation)
%          (optional argument, default=0.01)
%
% Outputs:
%    RH - Routh-Hurwitz matrix of size (length(P), length(P)/2).
%
% Example:
%     RH = rh([1, 2, 3, 4], 0.001);
%   Calculates the RH matrix for the polynomial s^3 + 2s^2 + 3s + 4. 
%   If one of the items on the first column from the matrix is zero, it will be replaced with the value eps=0.001.
%
% Author: Mircea Susca
% Technical University, Cluj-Napoca
% Email: mircea.m.susca@gmail.com
% March 2015; Last revision: 10-April-2015

%% parse eps
numvarargs = length(varargin);
if numvarargs > 1
    error('rh:TooManyInputs', ...
        'At most one extra argument.');
end

if isempty(varargin)
    eps = 0.01;
else
    eps = varargin{1};
end

%% Iterate through P and add alternatively on the first two rows of RH
len=length(P);
n=round(len/2);
RH=zeros(len,n);

i = 1;
j = 1;
for k = 1:len
    if rem(k,2)==1
        RH(1,j)=P(k);
        j=j+1;
    else
        RH(2,i)=P(k);
        i=i+1;
    end
end

%% Fill in the RH matrix
% - compute line values
% - null row => add coefficients from the above row's derivative polynomial 
% - first item null => replace with eps
for i=2:len
    prim=RH(i-1,1);
    for j=1:n-1
        if i == 2  % second row is not computed with the determinants
            break;
        end
        RH(i,j)=((RH(i-1,1)*RH(i-2,j+1))-(RH(i-2,1)*RH(i-1,j+1)))/prim;
    end
    
    if RH(i,:)==0
        grad=(len-i+1);
        dec=0;
        poz=1;
        for j=1:n-1
            RH(i,j)=(grad-dec)*(RH(i-1,poz));
            poz=poz+1;
            dec=dec+2;
        end
    end
    
    if RH(i,1)==0
        RH(i,1)=eps;
    end
end

%% Number of sign changes
p_rhp=0;	% number of right-half plane (RHP) poles
for i=2:len
    if RH(i-1,1)*RH(i,1) < 0    % <=> if sign(RH(i-1,1))*sign(RH(i,1)) == -1
        p_rhp = p_rhp + 1;
    end
end

fprintf('\n Number of right-half plane poles is %.0f.\n ',p_rhp);
if p_rhp == 0
    fprintf(' => System is stable.\n');
    is_stable = 1;
else
    fprintf(' => System is unstable.\n');
    is_stable = 0;
end

end