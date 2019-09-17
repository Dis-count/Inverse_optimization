function Inv_UFL3(fixed,v_i,u_ij,FC,TC)    % Inv_UFL1 中有四重循环，需要进一步优化 
% 该函数 功能 给定任意 m,n 可以得到UFL的 逆优化结果
% 解决了需要 n 重循环的问题
%  x_0 为给定服务对象以及是否开启 Facility  v_0 为原问题最优值  fixed 为目标值

% fixed = 29;   % 给定最优值  原问题最优值为29，因此 fixed = 29 时，得到最优结果是 0
 
% v_i = [0;0;1;1];
% 给定一个可行解
% u_ij = [0;0;0;0;
%         0;0;0;0;
%         1;0;1;0;
%         0;1;0;1;];

x_0 = [v_i;u_ij];

x0 =[x_0;x_0*(-1)];   % double the matrix    

% define primitive data
m = length(v_i);  
n = length(u_ij)/m; 

% 给出原问题的C_0

% FC  = [10; 10; 10; 10];

% M = 30;
% TC  = [ 3; 3; M; 2;
%          M; 1; 4; M;
%          3; M; 3; 4;
%          M; 2; M; 1;];
Costs =[FC;TC];    % Cost 为给定

V_0 = x_0'*Costs;    %  原最优解 


% Index helper function
% flowidx = @(w, p) nPlants * w + p;

% Build model
model.modelname = 'Inv_facility';
model.modelsense = 'min';

% Set data for variables
ncol = (m + m * n)*2 ;
model.lb    = zeros(ncol, 1);
model.ub    = [inf(ncol, 1)];
obj = [ones(m+m*n,1)];  

model.obj = [obj; obj];   % norm-1 c-Costs  均为正

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nplayers, 1)];

% Set data for constraints and matrix
nrow = m^n + 1;
model.A     = sparse(nrow, ncol);

model.sense = [repmat('=', 1, 1); repmat('>', m^n, 1)];

% Production constraints   注意限制条件需要遍历  这一点非常复杂


model.A(1,:) = x0;

model.rhs(1) = fixed- V_0;

a = eye(m);

% f = zeros(m^n,m+m*n);

s = 2;  


% varRange中数组取值对顺序敏感。
% 循环变量遍历集若是整体(而不是标量数字，或者单个字符等）构成，那么需要用花括号括起来。

% varRange = {1:4, 1:4, 1:4 , 1:4}; % 为按顺序的循环变量集合

b = repmat(1:m,1,n);   % n为player  即循环次数

varRange = mat2cell(b, [1], [repmat(m,1,n)]);

nloop = length(varRange);   % for循环层数，亦即循环变量个数

var = cell(nloop,1);    % 初始化循环变量

maxIter = cellfun(@length,varRange(:));

subs = arrayfun(@(x)1:x,maxIter,'un',0);

iter = cell(nloop,1);

[iter{end:-1:1}] = ndgrid(subs{end:-1:1});

temp = cellfun(@(x)x(:),iter,'un',0);

iter = num2cell(cat(2,temp{:}));    % 循环变量索引 

for k = 1:size(iter,1)  % 一次循环替代嵌套
    
    var = cellfun(@(x,y)x(y),varRange,iter(k,:),'un',0);
    
    sca = cell2mat(var);
    
    c = [a(:, sca)];   % 利用cell得到索引，从而根据索引取出独立的列向量
                    
    f = reshape(c',m*m,1); % 注意这里reshape是按列  而我们需要的变量顺序是按行因为将c转置 生成列向量
                
    % d = all(x<0)
    d = zeros(m,1);
    
    for j = 1:m            
        d(j) = 1-all(c(j,:)==0);   % 行中全为零则返回0
    end
    
    e = [d;f];   % 列向量
        
    ee = [e;e*(-1)];  % 必须为列向量
        
    model.A(s,:) = ee;
                    
    model.rhs(s) = fixed - e'* Costs;
        
    % f(s,:) = e;        
        
    s = s+1;
    
end


gurobi_write(model,'Inv_UFL3.lp');

% Guess at the starting point: close the plant with the highest fixed
% costs; open all others first open all plants
% model.start = [ones(nPlants, 1); inf(nPlants * nplayers, 1)];
% [~, idx] = max(FC);
% model.start(idx) = 0;



% Optimize
res = gurobi(model);

end