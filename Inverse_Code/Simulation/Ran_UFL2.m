% 本函数用于生成任意稀疏度 的 UFL逆优化 计算结果

% 给定任意稀疏度 如何生成 trans-cost 矩阵（即必须满足每一列都至少有一个不为M的值）是一个有意思的问题

% 本程序将之前所有的程序写在一起 成为 一个程序

%  本程序思路是 对于 每一类限制条件 使用一个向量来表示  

for m = 10:10:70
    myfun(m,2*m,50);
end

% This function is used to get excel to show the specific data result.
function myfun(m,n,k)   
% (m,n) 问题规模  k 为调用函数次数 默认 k=50;

theresult=zeros(k+1,3);

for j = 2:10   % j is multiplier

    for i=1:k 
    
        [a,b] = fun1(m,n,j);  % 重复调用  b 为最优值
        theresult(i,:) = [a,b,(a-0)/b];  % gap  最优  比值
         
    end
theresult(k+1,:) = [mean(theresult(1:k,1)),mean(theresult(1:k,2)),mean(theresult(1:k,3))];

filename = ['F:\Program Files\Matlab files\',num2str(m),'by',num2str(n),'-',num2str(j),'.xlsx'];

xlswrite(filename,theresult);

end

end




function [gap,opt1] = fun1(m,n,mul)  % 返回原UFL的最优值 和 本 LP 计算得到的 gap
% 给出 Facility 数量 m  Player 数量 n 

% mul 为 设施成本倍数

% 对应限制条件数量 为 （2mn+2m+2）个  变量有（5mn+3m+n) 个

fi = round(rand(1,m)*10*mul)';    % 可以调整 fi 与 rik 之间的比例 查看gap变化 
% rik = round(rand(1,m*n)*10)'; % 注意是列向量  （随机1-10,1-100,1-1000）
rik = fun3(m,n,round(m*n*0.25));


[opt1,opt2] = fun2(fi,rik);  % V_UFL 为 UFL最优值 
V_UFL = opt1;

% fi  = [10; 10; 10; 10];  
% 随机给出 facility cost 

% M = 30;    % define M= a bigger integer.


% 随机给出 transportation cost
% rik  = [ 3; 3; M; 2;
%          M; 1; 4; M;
%          3; M; 3; 4;
%          M; 2; M; 1;];


vi = opt2(1:m); 

uik = opt2(m+1:end);

% 产生一个 长度为 （5mn+3m+n) 的向量，每一个限制条件 产生 一个向量。    

% Build model
model.modelname = 'UFL';
model.modelsense = 'min';

% Set data for variables
ncol = 5*m*n + 3*m + n;

% 先试试变量大于零的情况 
model.lb    = zeros(ncol, 1);
model.ub    = inf(ncol, 1);
model.obj   = [zeros(n+m+3*m*n,1); ones(2*m + 2*m*n,1); ];

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
% gurobi_write(model,'UFL.lp');

% Optimize 
res = gurobi(model);

gap = res.objval;

end


function r = fun3(m,n,t)  % t为稀疏度
    x = ceil(rand(1,t)*m);    % 生成 t 个 1-m 的整数  作为坐标的行向量
    yy = 1:n;
    
    M = 100;     
    
    for i = 1:20
        y = ceil(rand(1,t)*n);  % 生成 t 个 1-n 的整数  作为坐标的列向量
        if sum(ismember(yy,y))== n  % 判断生成向量是否符合标准
        
            break
        else
            disp(i);    
        end
    end
    
    v = ceil(rand(1,t)*10);   % 生成 t 个 随机 坐标值
        
    r = full(sparse(x,y,v,m,n));
    s = sum(sum(r~=0));   % 判断 非零数 是否满足稀疏度
    if  s < t
        a = find(r==0); % 找到元素为 0 的位置
            
        for i = 1:t-s
            b = unidrnd(length(a));     % 随机取出 a 中一个数值
            r(a(b)) = ceil(rand(1,1)*10);    
            a(b) =[];   % 删除已添加的元素  注意这里 a 长度已经少 1 了
        end
    end
    r(r(:)==0) = M;    % 将稀疏矩阵中为 0值的地方 全设为 M
    r = reshape(r,m*n,1);
end


function [opt1,opt2] = fun2(fi,rik)   % fi rik 为列向量  求出 UFL的最优值以及向量
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


% Build model
model.modelname = 'I_UFL';
model.modelsense = 'min';

% Set data for variables
ncol = m*n + m ;

model.vtype = 'B';

model.obj   = [fi; rik];


% Set data for constraints and matrix

nrow = m*n + n;

model.A     = sparse(nrow, ncol);

model.rhs   = [zeros(m*n, 1); ones(n,1)];

model.sense = [repmat('>', m*n , 1); repmat('=', n , 1)];


for w =1:m
    for p =1:n
        model.A(p+n*(w-1),w) = 1;         % 第一组约束
        
        model.A(p+n*(w-1),m+p+(w-1)*n) = -1;
    end
end

for p = 1:n
    model.A(p+m*n, n*(0:(m-1))+m+p) = 1; % 第二组约束
end

% Save model
% gurobi_write(model,'I_UFL.lp');

% Optimize
% res = gurobi(model, params);
res = gurobi(model);
opt1 = res.objval;
opt2 = res.x;
% Print solution

end