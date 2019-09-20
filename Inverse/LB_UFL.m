% 这里有四组约束
% 第一个约束确保行内 r_ik (为1的最大值<为0 的最小值 )的大小关系

V_UFL = 28;   % 给定的最优值（目标值）

vi = [1 ; 1 ; 0;];

uik = [ 1; 0 ; 1 ;
        0; 1 ; 0 ;
        0; 0 ; 0 ;];

m = length(vi);
n = length(uik)/m;

x_0 = [vi; uik];
x0 =[x_0; x_0*(-1)];  %求解需要的double向量

% 给定的原问题Costs

FC    = [5; 6; 7;];

TC    = [11; 4; 8;
         5; 7; 10;
         19; 6; 3;];

Costs =[FC;TC];

V_0 = x_0'*Costs;    %  给定的最优解向量 对应的原成本矩阵的花费值

v1 = find(vi == 1);   % v1记录 vi=1 的横坐标 的 列向量   要分别找到vi为 0 或 1 的下标

v0 = find(vi == 0);

v1num = length(v1);

u1 = cell(v1num,3);  % 只记录 不确定长度的 vi=1 的 位置及相应的TC值 需要删掉剩下 v1num-sss

u11 = zeros(v1num,5);  % 记录最大最小坐标

ui = reshape(uik,m,n)';   % 转置 得到uik 的 解矩阵

TCi = reshape(TC,m,n)';   % 得到矩阵

s = 0;  % 用于记录 vi 中 rik==1 个数为1 的数量

sss = 0;

v1min = max(TCi(1,:));  % 记录 全局最小值

Loc =zeros(1,2);  % 记录 最小值位置  约束2用不到了

loc = zeros(v1num,3);  % 记录只有一个vi的 行 纵 坐标 和 值  需要删掉剩下 v1num-s 个
% 约束3用到   约束2也用到

t=1;

for i = v1'  % v1 需要是行向量 % 循环 每一行

    a1 = find(ui(i,:) == 1);
    a0 = find(ui(i,:) == 0);

    b = TCi(i,:);  % 找出TC中对应 ui 为 1 的部分

    [umax,in] = max(b(a1));

    [umin,ni] = min(b(a0));

    % a1(in)  得到筛选过后 的 最大值坐标索引

    u11(t,:) = [i, umax, a1(in), umin, a0(ni)];  % 里面存放 vi为 1的 坐标 对应行的 最小值和最大值

    if length(a1)==1

        s = s + 1;

        loc(s,:) = [i,a1,b(a1)];   % 记录只有一个 vi=1 的 行纵值

    else
        sss = sss + 1;  % 1的值为两个及以上

        u1{sss,1} = i;  % 记录 行号

        u1{sss,3} = b(a1);  % 第三列 存储 uik 为 1 的部分对应 TC 的值

        u1{sss,2} = a1;  % 单元数组 第二列里面存放中 rik 对应行中 为1 的坐标

    end

    if umin < v1min

      Loc = [i,a0(ni)];  % 记录  行 和 列

      v1min = umin;  % 更新全局 uik 最小值

    end

    t=t+1;
end

% num_row =  cell2mat(u1(:,1));   % 记录行号 形成向量

% num_col =  cell2mat(u1(:,2));   % 记录列值 形成向量

% num_value = cell2mat(u1(:,3));  % 记录 相应 TC

% 删除 0 向量
loc(s+1:end,:)=[];
u1(sss+1:end,:)=[];

u0 = zeros(m-v1num,5);  % 此处 重新定义 一个类似的单元数组 存放为0 的部分

t = v1num + 1;

for i = v0'  % v0 需要是行向量

    [umax,in] = max(TCi(i,:));

    [umin,ni] = min(TCi(i,:));

    u0(t-v1num,:) = [i,umax, in, umin, ni];  % 单元数组 第一列里面存放 vi为 1的坐标

    t=t+1;
end

% 需要给出 vi 为1的下标  以及 为0 的下标

% [h l] = max(reshape(TC,m,n), [], 2);   % 给出每一行的最大值列向量h 以及 下标向量l。  注意1 为每列， 2 为每行。

model.modelname = 'LB_Inv_UFL';
model.modelsense = 'min';

col = m + m * n;
ncol = col * 2 ;

model.lb  = zeros(ncol, 1);
model.ub  = inf(ncol, 1);

obj = ones(col ,1);

model.obj = [obj; obj];   % norm-1 c-Costs  均为正

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nplayers, 1)];

% Set data for constraints and matrix

% nrow = m + n + v1num - s + 1; % 前两个约束加起来为 m 个；第三个约束有 n 个；
% 第四个约束有 vi=1 对应行中 rik==1 个数大于 1 的个数 加 一个等式

nrow = v1num + (m - v1num)*n + n + 1 + (v1num-s) ;

model.A     = sparse(nrow, ncol);

model.sense = [repmat('=', 1, 1); repmat('<', nrow-1, 1)];

% Production constraints   注意限制条件需要遍历  这一点非常复杂

model.A(1,:) = x0;

model.rhs(1) = V_UFL-V_0;


%  第一类约束
for i = 1 : v1num
    % 即需要找到 i 所在行的 uik 为 1 的最大值  % 需要给出下标

    inde = zeros(col,1);

    % u1{t,1}()  从 cell 中取出元素

    inde((u11(i,1)-1)*n + m + u11(i,3)) = 1;  % max-index

    inde((u11(i,1)-1)*n + m + u11(i,5)) = -1; % min-index

    indexd = [inde; inde*(-1)];

    rhs = (-1)*u11(i,2) + u11(i,4);

    model.A(1 + i, :) = indexd;  % 有问题

    model.rhs(1 + i) = rhs;

end

%  第二类约束   % 这里有  (m-v1num)*s个
for k= 1:(m - v1num)

  for i=1:s

    inde = zeros(col,1);

    inde(u0(k,1)) = -1; % fi==0的横坐标

    a = v0(k);

    inde((u0(k,1)-1)*n + m + loc(i,2)) = -1;  % 使用 s 的 纵坐标

    inde((loc(i,1)-1)*n + m + loc(i,2)) = 1;  % 对应列 rik

    inde(loc(i,1)) = 1;

    indexd = [inde; inde*(-1)];

    rhs = TC((u0(k,1)-1)*n + loc(i,2)) + FC(a) - FC(loc(i,1)) - loc(i,3);

    model.A(v1num + 1 + i + (k-1)*s, :) = indexd;

    model.rhs(v1num + 1 + i + (k-1)*s) = rhs;

  end

end


for k= 1:(m-v1num)  % 这里有  (m-v1num)*(n-s) 个

  tt = 1;

  for i = 1 : sss  % 行数

    t=1;

    for j = cell2mat(u1(i,2))

      inde = zeros(col,1);

      inde(u0(k,1)) = -1; % fi==0的横坐标

      a = v0(k);

      inde((u0(k,1)-1)*n + m + j) = -1;  % 使用 s 的 纵坐标

      inde((cell2mat(u1(i,1))-1)*n + m + j) = 1;  % 对应列 rik

      indexd = [inde; inde*(-1)];

      num_value = cell2mat(u1(i,3));

      rhs = TC((u0(k,1)-1)*n + j) + FC(a) - num_value(t);


      model.A(v1num + (m-v1num)*s + 1 + tt + (k-1)*(n-s), :) = indexd;

      model.rhs(v1num + (m-v1num)*s + 1 + tt + (k-1)*(n-s)) = rhs;

      t = t + 1;

      tt = tt + 1;

    end

  end

end

% 第三类约束   % 这里有  (s) 个
for i = 1 : s

   inde = zeros(col,1);

   inde(loc(i,1)) = 1;   % 有问题

   inde((loc(i,1)-1)*n + m + loc(i,2)) = 1;

   mmi = cell2mat(u1(:,1));

   [mi,mindex] = min(TCi(mmi,loc(i,2))); % 最小值索引即为行号

   inde((mmi(mindex)-1)*n + m + loc(i,2)) = -1;  %index 有问题

   indexd = [inde; inde*(-1)];

   rhs = -loc(i,3) - FC(loc(i,1)) + mi;

   model.A(v1num + (m-v1num) * n + 1 + i, :) = indexd;

   model.rhs(v1num + (m-v1num) * n + 1 + i) = rhs;

end

tt = 1;

for i = 1 : sss  % sss行数   % 这里有 (n-s) 个

  t = 1;

  for j = cell2mat(u1(i,2))

    inde = zeros(col,1);

    inde((cell2mat(u1(i,1))-1)*n + m + j) = 1;   % 对应列 rik

    num_value = cell2mat(u1(i,3));

    mmi = loc(:,1);

    [mi,mindex] = min(TCi(mmi, j));   % 一个的行  多个的对应列

    inde((mmi(mindex)-1)*n + m + j) = -1;

    rhs = mi - num_value(t) ;

    indexd = [inde; inde*(-1)];

    model.A(v1num + (m-v1num) * n + 1 + s + tt, :) = indexd;

    model.rhs(v1num + (m-v1num) * n + 1 + s + tt) = rhs;

    t = t + 1;

    tt = tt + 1;

  end

end

% 第四类约束
for i = 1 : (v1num-s)

  inde = zeros(col,1);

  inde(cell2mat(u1(i,1))) = 1;

  in = cell2mat(u1(i,2))+ m + (cell2mat(u1(i,1))-1)*n; %  位置

  % cell2mat(u1(i,2)) 为 列 形成的向量

  inde(in) = 1;

  new = TCi(loc(1:s,1),cell2mat(u1(i,2)));   % 找到 相应列

  aa = loc(:,1);

  [mii,ind] = min(sum(new,2));  % 对 行 求和  得到 最小值 和 行索引

  inde((aa(ind)-1)*n + m + cell2mat(u1(i,2))) = -1;

  indexd = [inde; inde*(-1)];

  rhs = - sum(cell2mat(u1(i,3))) - FC(cell2mat(u1(i,1))) + mii ;

  model.A(v1num + (m-v1num) * n + 1 + n + i, :) = indexd;

  model.rhs(v1num + (m-v1num) * n + 1 + n + i) = rhs;

end

gurobi_write(model,'LB_UFL.lp');

res = gurobi(model);
