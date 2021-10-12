function [x,status] = LPSol(A,b,c)

% code minimizes c^T x 
  % subject to: x >=0
%               A*x = b
%
% this function solves an LP
%
[m,n]    = size(A);
msk      = find(b < 0);
b(msk)   = - b(msk);
A(msk,:) = - A(msk,:);
msk      = find(b == 0);
if (length(msk)>0)
    disp(['LP is degenerate']);
    status.feas = 'degenerate';
    x           = NaN*ones(n,1);
    return;
end
% 
% Phase I
%
Ahat     = [A, eye(m)];
chat     = zeros(m+n,1);
chat(n+1:n+m)  = ones(m,1);
Bhat           = transpose((n+1:n+m));
[xhat,yhat,B,simplex]  = Simplex(Ahat,b,chat,Bhat);
if (simplex == 1)
    disp(['LP is degenerate']);
    status.feas = 'degenerate';
    x           = NaN*ones(n,1);
    return;
end
disp('simplex method successful for Phase I')
if (max(B)>n)
    disp(['LP is infeasible, Phase I objective = ', num2str(transpose(chat)*xhat)])
    status.feas = 'infeasible';
    x = [];
    return;
end
% 
% Phase II
%
[x,y,B,simplex]  = Simplex(A,b,c,B);
if (simplex == 1)  
    disp(['LP is feasible but degenerate']);
    status.feas = 'degenerate';
    x           = NaN*ones(n,1);
    return;
end
if (simplex == -1)
    disp(['LP is feasible but without optimal']);
    status.feas = 'unbounded';
    status.B    = B;
    return;
end
disp('simplex method successful for Phase II')
disp(['LP optimal reached at objective = ', num2str(transpose(c)*x)])
status.feas = 'optimal';
status.B    = B;

