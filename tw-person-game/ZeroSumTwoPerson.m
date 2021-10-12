%
% This code solves the zero sum, two person game 
%

function   [p,q,val] = ZeroSumTwoPerson(A)
%
%
[m,n] = size(A);
em    = ones(m,1);
en    = ones(1,n);
%
% put game in canonical LP form
%
AA    = [A, -em, em, eye(m); en, zeros(1,2+m)];
bb    = zeros(m+1,1);
bb(m+1) = 1;
bb(1:m) = 1e-6 * abs(randn(m,1));
cc      = zeros(m+n+2,1);
cc(n+1) =  1;
cc(n+2) = -1;

%
% Simplex Method
%
[x,status] = LPSol(AA,bb,cc);
feas       = status.feas;
opt        = 'optimal';
success    = 1;
if (length(feas)==length(opt))
    if (feas == opt)
        success    = 0;
    end
end
if (success >0)
    disp('solver failure')
    q = [];
    p = [];
    val = NaN;
else
    disp('optimal strategy found')
    bb(1:m)   = 0;
    B         = status.B;
    x(B)      = (AA(:,B)) \ bb;
    q         = x(1:n);
    y         = (transpose(AA(:,B))) \cc(B);
    p         =-y(1:m);
    p(p < 1e-15) = 0;
    q(q < 1e-15) = 0;
    val       = y(m+1);
end
