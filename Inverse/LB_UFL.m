% 这里有四组约束
% 第一个约束确保行内 r_ik (为1的最大值<为0 的最小值 )的大小关系

V_UFL = 28;   % 给定的最优值（目标值）

vi = [1 ; 1 ; 0;];

uik = [ 1; 0 ; 1 ;
        0; 1 ; 0 ;
        0; 0 ; 0 ;];

m = length(vi);
n = length(uik)/m;

x_0 = [vi;uik];
x0 =[x_0;x_0*(-1)];  %求解需要的double向量

% 给定的原问题Costs

FC    = [5; 6; 7;];

TC    = [11; 4; 8;
         5; 7; 10;
         19; 6; 3;];

Costs =[FC;TC];

V_0 = x_0'*Costs;    %  给定的最优解向量 对应的原成本矩阵的花费值

v1 = find(vi == 1);   % v1记录 vi=1 的横坐标 的 列向量   要分别找到vi为 0 或 1 的下标

v0 = find(vi == 0);

u1 = cell(length(v1),2);

ui = reshape(uik,m,n)';   % 转置 得到uik 的 解矩阵

TCi = reshape(TC,m,n)';   % 得到矩阵

s = 0; %  用于记录 vi 中 rik==1 个数为1 的数量

for i = v1'  % v1 需要是行向量

    t=1;

    [umax,in] = max(TCi(i,:));

    [umin,ni] = min(TCi(i,:));

    u1{t,1} = [i,[umax, in, umin, ni]];  % 单元数组 第一列里面存放 vi为 1的 坐标 对应行的 最小值和最大值

    u1{t,2} = find(ui(i,:) == 1);  % 单元数组 第二列里面存放中 rik 对应行中 为1 的坐标

    if length(find(ui(i,:) == 1))==1

        s=s+1;

    end

    t=t+1;
end

u0 = cell(m-length(v1),2);  % 此处 重新定义 一个类似的单元数组 存放为0 的部分

for i = v0'  % v1 需要是行向量
    t=1;

    [umax,in] = max(TCi(i,:));

    [umin,ni] = min(TCi(i,:));

    u0{t,1} = [i,[umax, in, umin, ni]];  % 单元数组 第一列里面存放 vi为 1的坐标
    u0{t,2} = find(ui(i,:) == 0);  %单元数组 第二列里面存放中 rik 对应行中 为0 的坐标

    t=t+1;
end

fi;

% 需要给出 vi 为1的下标  以及 为0 的下标

rik;

[h l] = max(reshape(TC,m,n), [], 2);   % 给出每一行的最大值列向量m 以及 下标向量l。  注意1 为每列， 2 为每行。

model.modelname = 'LB_Inv_UFL';
model.modelsense = 'min';

col = m + m * n;
ncol = col * 2 ;

model.lb    = zeros(ncol, 1);
model.ub    = inf(ncol, 1);

obj = ones(col ,1);

model.obj = [obj; obj];   % norm-1 c-Costs  均为正

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nplayers, 1)];

% Set data for constraints and matrix

nrow = m + n + length(v1) - s + 1  ; % 前两个约束加起来为 m 个；第三个约束有 n 个；
% 第四个约束有 vi=1 对应行中 rik==1 个数大于 1 的个数 加 一个等式

model.A     = sparse(nrow, ncol);

model.sense = [repmat('=', 1, 1); repmat('<', nrow-1, 1)];

% Production constraints   注意限制条件需要遍历  这一点非常复杂

model.A(1,:) = x0;

model.rhs(1) = V_UFL-V_0;


%  第一个约束
for i=1:length(v1)
    % 即需要找到 i 所在行的 uik 为 1 的最大值  % 需要给出下标

    index = zeros(col,1);

    % u1{t,1}()  从 cell 中取出元素

    index((u1{i,1}(1)-1)*n + m + u1{i,1}(3)) = 1;  % max-index
    index((u1{i,1}(1)-1)*n + m + u1{i,1}(5)) = -1; % min-index

    indexd = [index;index*(-1)];

    rhs = (-1)* u1{i,1}(2)+ u1{i,1}(4);

    model.A(1+i,:) = indexd;

    model.rhs(1+i) = rhs;

end

%  第二个约束
for i=1:length(v0)

    index = zeros(col,1);

    index(v0(i)) = -1 % fi==0
    index((u1{i,1}(1)-1)*n + m + u1{i,1}(5)) = -1


    index((u1{i,1}(1)-1)*n + m + u1{i,1}(3)) = 1;

    indexd = [index;index*(-1)];

    rhs = (-1)* u1{i,1}(2)+ u1{i,1}(4);

    model.A(1+i,:) = indexd;

    model.rhs(1+i) = rhs;


end
