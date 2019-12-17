function [opt1,opt2] = CFL(FC,TC,k,d)
% FC TC 为列向量  求出给定 costs ��?  UFL的最优�??
% 给出 Facility 数量 m  Player 数量 n

% The basic CFL Optimization solution.

% 对应限制条件数量 ��? （mn+m）个  变量有（mn+m) ��?

m = length(FC) ;

n = length(TC)/m ;

% FC  = [10; 10; 10; 10];

% M = 100;    % define M= a bigger integer.

% 注意这里应该使用循环生成向量��? 但这里为了计算简单的例子（n=4)，我们直接手动添加变量�??

% TC  = [ 3; 3; M; 2;
%          M; 1; 4; M;
%          3; M; 3; 4;
%          M; 2; M; 1;];

% 产生��?��? 长度��? （mn+m) 的向量，每一个限制条��? 产生 ��?个向量�??

% Build model
model.modelname = 'CFL';
model.modelsense = 'min';

% Set data for variables
ncol = m*n + m ;

% model.vtype = 'B';

model.obj   = [FC; TC];

model.vtype = [repmat('B', m, 1); repmat('C', m*n, 1)];
%
% for p = 1:nPlants
%     model.varnames{p} = sprintf('Open%d', p);
% end
%
% for w = 1:nWarehouses
%     for p = 1:nPlants
%         v = flowidx(w, p);
%         model.varnames{v} = sprintf('Trans%d,%d', w, p);
%     end
% end

% Set data for constraints and matrix

nrow = m + n;

model.A     = sparse(nrow, ncol);

model.rhs   = [zeros(m, 1); ones(n,1)];

model.sense = [repmat('>', m , 1); repmat('=', n , 1)];

for w =1:m                       % 第一组约��?

    model.A(w,w) = k(w);

    for p =1:n

        model.A(w, m+(w-1)*n+p) = -d(p);
        %��?定要注意这里��? m+m*n��?  注意循环顺序
    end

end


for p = 1:n

    model.A(p+m, n*(0:(m-1))+m+p) = 1; % 第二组约��?

end

% Save model
% gurobi_write(model,'UFL.lp');

% Guess at the starting point: close the plant with the highest fixed
% costs; open all others first open all plants

% model.start = [ones(m, 1); inf(nPlants * nWarehouses, 1)];
% [~, idx] = max(FixedCosts);
% model.start(idx) = 0;

% Set parameters
% params.method = 2;

% Optimize
% res = gurobi(model, params);

params.outputflag = 0;

res = gurobi(model,params);

opt1 = res.objval;   % Solution

opt2 = res.x;        % Value

% Print solution

end
