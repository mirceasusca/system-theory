function [is_stable,c1]=jury_criterion(p)
%% 
% Apply the Schur-Cohn-Jury stability test on characteristic polynomial p:
% p(z) = a(0)*z^n + a(1)*z^(n-1) + ... + a(n-1)*z + a(n)
%    
%    col.1  col.2  col.3  ...  col.n-1 col.n
%   -------------------------------------------
%  |  a(0)   a(1)   a(2)   ...  a(n-1)  a(n)
%  |  a(n)  a(n-1) a(n-2)  ...   a(1)   a(0)
%  |-------------------------------------------
%  |  b(0)   b(1)   b(2)   ...  b(n-1)
%  | b(n-1) b(n-2) b(n-3)  ...   b(0)
%  |-------------------------------------------
%  |  c(0)   c(1)   c(2)   ...
%  | c(n-2) c(n-3) c(n-3)  ... 
%   -------------------------------------------
%  
%  where b(k) = a(k) - a(n)/a(0) * a(n-k) etc.
% 
% If all the terms in the first column of the odd rows are positive, the
% polynomial p(z) is stable.
%
% Example:
%   h=zpk([],[0.98, 0.235, -0.923],1,1);
%   [num,den]=tfdata(h,'v');
%   [is_stable,c1]=jury_criterion(denz);
%   >> is_stable = 1; c1 = [0.9548,0.2243,0.0168]
%
% Author: Mircea Susca
% Technical University, Cluj-Napoca
% Email: mircea.m.susca@gmail.com
% 10-Dec-2017.
%

%%
if p(1) == 0
    error('Error. The first element of the array must be non-zero');
end

%% Initialization
p = p/sign(p(1)); % normalize polynomial to positive dominant term
c1 = zeros(1,length(p)-1);  % first column of the Jury table

r = p; % the first odd row starts with the polyonomial p(z)

%%
for i=1:length(p)-1
    n = length(r);
    rp = zeros(1,n-1);  % initialization for the next odd row
   
    % fill in the next row
    for j=1:n-1
        rp(j) = r(1)*r(j) - r(n)*r(n+1-j);
    end
   
    c1(i) = rp(1);  % add the element to the column
    r = rp;
end

if min(c1) > 0
    is_stable = true;
else
    is_stable = false;
end