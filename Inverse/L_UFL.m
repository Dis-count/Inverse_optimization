%  本程序思路是 对于 每一类限制条件 使用一个向量来表示  
function L_UFL(m,n)
% 给出 Facility 数量 m  Player 数量 n 

% 对应限制条件数量 为 （2mn+2m+2）个  变量有（5mn+3m+n) 个

V_UFL = 28;  % V_UFL 为 给定 目标值


%  （fi,rik)为 原向量 C_0
% fi  = [10; 10; 10; 10];   
fi  = [5; 6; 7;];


% M = 30;    % define M= a bigger integer.

% 注意这里应该使用循环生成向量； 但这里为了计算简单的例子（n=4)，我们直接手动添加变量。

% rik  = [ 3; 3; M; 2;
%          M; 1; 4; M;
%          3; M; 3; 4;
%          M; 2; M; 1;];
rik  = [11; 4; 8; 
         5; 7; 10; 
         19; 6; 3; ];

vi = [1 ; 1 ; 0;];

uik = [ 1; 0 ; 1 ; 
        0; 1 ; 0 ;  
        0; 0 ; 0 ;];
     
    
     
% vi = [0 ; 0 ; 1 ; 1]; 
% 
% uik = [ 0 ; 0 ; 0 ; 0 ;
%         0 ; 0 ; 0 ; 0 ;
%         1 ; 0 ; 1 ; 0 ; 
%         0 ; 1 ; 0 ; 1 ;];
% 这里给出的是原优化问题的最优解 但只要是一个可行解就可以。
% 产生一个 长度为 （6mn+4m+n) 的向量，每一个限制条件 产生 一个向量。
    

% Facility location: a company currently ships its product from 5 plants
% to 4 warehouses. It is considering closing some plants to reduce
% costs. What plant(s) should the company close, in order to minimize
% transportation and fixed costs?
%
% Note that this example uses lists instead of dictionaries.  Since
% it does not work with sparse data, lists are a reasonable option.


% define primitive data
% nPlants     = m;
% nWarehouses = n;


% Fixed costs for each plant
% FixedCosts  = [12000; 15000; 17000; 13000; 16000];

% Transportation costs per thousand units
% TransCosts  = [
%     4000; 2000; 3000; 2500; 4500;
%     2500; 2600; 3400; 3000; 4000;
%     1200; 1800; 2600; 4100; 3000;
%     2200; 2600; 3100; 3700; 3200];

% Index helper function
%flowidx = @(w, p)  n* p + w;

% Build model
model.modelname = 'UFL';
model.modelsense = 'min';

% Set data for variables
ncol = 5*m*n + 3*m + n;

% 先试试变量大于零的情况 
model.lb    = zeros(ncol, 1);
model.ub    = [inf(ncol, 1)];
model.obj   = [zeros(n+m+3*m*n,1); ones(2*m + 2*m*n,1); ];
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
nrow = 2*m*n+2*m+2;

model.A     = sparse(nrow, ncol);


model.rhs   = [V_UFL; zeros(m + m*n, 1); V_UFL; fi; rik];
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
gurobi_write(model,'UFL.lp');


% Guess at the starting point: close the plant with the highest fixed
% costs; open all others first open all plants

% model.start = [ones(m, 1); inf(nPlants * nWarehouses, 1)];
% [~, idx] = max(FixedCosts);
% model.start(idx) = 0;

% Set parameters
% params.method = 2;

% Optimize
% res = gurobi(model, params);
gurobi(model);

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

