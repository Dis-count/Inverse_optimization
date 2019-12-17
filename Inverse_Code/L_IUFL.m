%  本程序思路是 对于 每一类限制条件 使用一个向量来表示
% 将具体例子抽象化
function res = L_IUFL(V_UFL,vi,uik,FC,TC)

% V_UFL = 28;  % V_UFL 为 给定 目标值

% FC  = [5; 6; 7;];

% 注意这里应该使用循环生成向量； 但这里为了计算简单的例子（n=4)，我们直接手动添加变量。

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
% 这里给出的是原优化问题的最优解 但只要是一个可行解就可以。
m = length(vi);

n = length(uik)/m;

% Build model
model.modelname = 'L_IUFL';
model.modelsense = 'min';

% Set data for variables
ncol = 5*m*n + 3*m + n;

% 先试试变量大于零的情况

model.lb    = zeros(ncol, 1);
model.ub    = [inf(ncol, 1)];
model.obj   = [zeros(n+m+3*m*n,1); ones(2*m + 2*m*n,1);];

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nWarehouses, 1)];
%

% Set data for constraints and matrix
nrow = 2*m*n+2*m+2;

model.A     = sparse(nrow, ncol);

model.rhs   = [V_UFL; zeros(m + m*n, 1); V_UFL; FC; TC];
model.sense = [repmat('>', 1, 1); repmat('=', 2*m*n + 2*m + 1, 1)];

model.A(1,1:n) = 1;  % 第一类约束

for p = 1:m
    for w = 1:n
        model.A(p+1, n*p+w) = 1;
    end
    model.A(p+1, n+n*m+p) = -1;
%    model.constrnames{p} = sprintf('Capacity%d', p);
end

% 第二类约束
for p = 1:m
    for w = 1:n
        model.A((p-1)*n+w+m+1,[w,m*n+m+p*n+w]) = 1;

        model.A((p-1)*n+w+m+1,[p*n+w,2*m*n+m+p*n+w]) = -1;
    end
end   % 第三类约束

for p = 1:m
    for w = 1:n
    model.A(m*n+m+2,m+2*m*n+p*n+w) = uik((p-1)*n+w);

    end
    model.A(m*n+m+2,n+m*n+p) = vi(p);
%    model.constrnames{nPlants+w} = sprintf('Demand%d', w);
end   % 第四个约束


for p = 1:m
    model.A(m*n+m+2+p,n+m+3*m*n+p) = -1;
    model.A(m*n+m+2+p,[n+m*n+p,n+2*m+3*m*n+p]) = 1;   % 保持右侧约束为正的fi 下同
end   % 第五个约束


for p = 1:m
    for w = 1:n
        model.A((p-1)*n+m*n+2*m+2+w, 3*m+3*m*n+p*n+w) = -1;

        model.A((p-1)*n+m*n+2*m+2+w, [m+2*m*n+p*n+w,3*m+4*m*n+p*n+w]) = 1;
    end
end  % 第六个约束

% Save model
% gurobi_write(model,'UFL.lp');

% Optimize
% res = gurobi(model, params);

result = gurobi(model);

res = result.objval;

end
