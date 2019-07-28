  
function [opt1,opt2] = I_UFL(fi,rik)   % fi rik 为列向量   求出给定 costs 时  UFL的最优值
% 给出 Facility 数量 m  Player 数量 n 

% 对应限制条件数量 为 （mn+m）个  变量有（mn+m) 个

m = length(fi) ;
n = length(rik)/m ;


% disp(fi');

%fprintf('Facility Costs: %g\n', fi);

% disp(rik');

%fprintf('Trans Costs: %f\n', rik);

% fi  = [10; 10; 10; 10];   

% M = 100;    % define M= a bigger integer.

% 注意这里应该使用循环生成向量； 但这里为了计算简单的例子（n=4)，我们直接手动添加变量。


% rik  = [ 3; 3; M; 2;
%          M; 1; 4; M;
%          3; M; 3; 4;
%          M; 2; M; 1;];


% 产生一个 长度为 （mn+m) 的向量，每一个限制条件 产生 一个向量。
    

% Facility location: a company currently ships its product from 5 plants
% to 4 warehouses. It is considering closing some plants to reduce
% costs. What plant(s) should the company close, in order to minimize
% transportation and fixed costs?
%
% Note that this example uses lists instead of dictionaries.  Since
% it does not work with sparse data, lists are a reasonable option.


% Index helper function
%flowidx = @(w, p)  n* p + w;

% Build model
model.modelname = 'I_UFL';
model.modelsense = 'min';

% Set data for variables
ncol = m*n + m ;

model.vtype = 'B';

model.obj   = [fi; rik];

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nWarehouses, 1)];
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

nrow = m*n + n;

model.A     = sparse(nrow, ncol);

model.rhs   = [zeros(m*n, 1); ones(n,1)];

model.sense = [repmat('>', m*n , 1); repmat('=', n , 1)];


for w =1:m
    for p =1:n
        model.A(p+n*(w-1),w) = 1;         % 第一组约束
        
        model.A(p+n*(w-1),n*w+p) = -1;
    end
end

for p = 1:n
    model.A(p+m*n, [n*(0:(m-1))+n+p]) = 1; % 第二组约束
end

% Save model
% gurobi_write(model,'I_UFL.lp');


% Guess at the starting point: close the plant with the highest fixed
% costs; open all others first open all plants

% model.start = [ones(m, 1); inf(nPlants * nWarehouses, 1)];
% [~, idx] = max(FixedCosts);
% model.start(idx) = 0;

% Set parameters
% params.method = 2;

% Optimize
% res = gurobi(model, params);
res = gurobi(model);
opt1 = res.objval;
opt2 = res.x;
% Print solution

% if strcmp(res.status, 'OPTIMAL')
%     fprintf('\nTotal Costs: %g\n', res.objval);
%     fprintf('solution:\n');
%     for p = 1:nPlants
%         if res.x(p) > 0.99
%             fprintf('Plant %d open:\n', p);
%         end
%         for w = 1:nWarehouses
%             if res.x(flowidx(w, p)) > 0.0001
%                 fprintf('  Transport %g units to warehouse %d\n', res.x(flowidx(w, p)), w);
%             end
%         end
%     end
% else
%     fprintf('\n No solution\n');
% end

end

