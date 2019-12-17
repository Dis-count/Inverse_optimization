%  This program is designed by every constraints is expressed by a vector.

function res = H_ICFL(V_UFL,vi,uik,FC,TC)
% This is a heuristic method to calculate the inverse CFL.
% V_UFL = 28;  % V_UFL 为 给定 目标值
% k_i is the capacity of the facility.
% d_i is the demand of the customer.
% FC  = [5; 6; 7;];

% TC  = [ 3; 3; M; 2;
%          M; 1; 4; M;
%          3; M; 3; 4;
%          M; 2; M; 1;];

% vi = [0 ; 0 ; 1 ; 1];
%
% uik = [ 0 ; 0 ; 0 ; 0 ;
%         0 ; 0 ; 0 ; 0 ;
%         1 ; 0 ; 1 ; 0 ;
%         0 ; 1 ; 0 ; 1 ;];
% For the convinience of calculation.
% 这里给出的是原优化问题的最优解 但只要是一个可行解就可以。
m = length(vi);

n = length(uik)/m;

k = ones(m,1);

d = ones(n,1);

% Build model
model.modelname = 'H_ICFL';
model.modelsense = 'min';

% Set the number of variables

ncol = 4*m*n + 4*m + n;

% variable > 0

model.lb    = zeros(ncol, 1);
model.ub    = [inf(ncol, 1)];
model.obj   = [zeros(n+2*m+2*m*n,1); ones(2*m + 2*m*n,1);];

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nWarehouses, 1)];
%

% Set the number of constraints.

nrow = 2*m*n+2*m+2;

model.A     = sparse(nrow, ncol);

model.rhs   = [V_UFL; zeros(m + m*n, 1); V_UFL; FC; TC];

model.sense = [repmat('>', 1, 1); repmat('=', 2*m*n + 2*m + 1, 1)];

model.A(1,1:n) = 1;  % The first class of constraints

for p = 1:m

        model.A(p+1, n+p) = k(i);

    end

    model.A(p+1, n+m+p) = -1;

end

% The Second class of constraints

for p = 1:m
    for w = 1:n
        model.A((p-1)*n+w+m+1,[w,2*m+p*n+w]) = 1;

        model.A((p-1)*n+w+m+1,m*n+2*m+p*n+w) = -1;

        model.A((p-1)*n+w+m+1,n+p) = -d(k);

    end
end   % The Third class of constraints


for p = 1:m
    for w = 1:n

    model.A(m*n+m+2,2*m+*m*n+p*n+w) = uik((p-1)*n+w);

    end

    model.A(m*n+m+2,n+m+p) = vi(p);

end   % The Fourth class of constraints


for p = 1:m
    model.A(m*n+m+2+p,n+2*m+2*m*n+p) = -1;

    model.A(m*n+m+2+p,[n+m+p,n+3*m+2*m*n+p]) = 1;   % 保持右侧约束为正的fi 下同
end   % The Fifth class of constraints


for p = 1:m
    for w = 1:n
        model.A((p-1)*n+m*n+2*m+2+w, 4*m+2*m*n+p*n+w) = -1;

        model.A((p-1)*n+m*n+2*m+2+w, [2*m+m*n+p*n+w,4*m+3*m*n+p*n+w]) = 1;
    end
end  % The Sixth class of constraints

% Save model
% gurobi_write(model,'UFL.lp');

% Optimize
% res = gurobi(model, params);

result = gurobi(model);

res = result.objval;

end
