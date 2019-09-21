% 这里有四组约��?
% 第一个约束确保行��? r_ik (��?1的最大�??<��?0 的最小�?? )的大小关��?

function res = LB_IUFL(V_UFL,vi,uik,FC,TC)

% V_UFL = 28;   % 给定的最优�?�（目标值）

% vi uik ��? 给定��?优解

% FC TC  ��? 原矩阵�??

m = length(vi);

n = length(uik)/m;

x_0 = [vi; uik];

x0 = [x_0; x_0*(-1)];  %求解��?要的double向量

% 给定的原问题Costs

Costs =[FC;TC];

V_0 = x_0'*Costs;    %  给定的最优解向量 对应的原成本矩阵的花费�??

v1 = find(vi == 1);   % v1记录 vi=1 的横坐标 ��? 列向��?   要分别找到vi��? 0 ��? 1 的下��?

v0 = find(vi == 0);

v1num = length(v1);

u1 = cell(v1num,3);  % 只记��? 不确定长度的 vi=1 ��? 位置及相应的TC��? ��?要删掉剩��? v1num-sss

u11 = zeros(v1num,5);  % 记录��?大最小坐��?

ui = reshape(uik,m,n)';   % 转置 得到uik ��? 解矩��?

TCi = reshape(TC,m,n)';   % 得到矩阵

s = 0;  % 用于记录 vi ��? rik==1 个数��?1 的数��?

sss = 0;

v1min = max(TCi(1,:));  % 记录 全局��?小�??

Loc =zeros(1,2);  % 记录 ��?小�?�位��?  约束2用不到了

loc = zeros(v1num,3);  % 记录只有��?个vi��? ��? ��? 坐标 ��? ��?  ��?要删掉剩��? v1num-s ��?
% 约束3用到   约束2也用��?

t=1;

for i = v1'  % v1 ��?要是行向��? % 循环 每一��?

    a1 = find(ui(i,:) == 1);
    a0 = find(ui(i,:) == 0);

    b = TCi(i,:);  % 找出TC中对��? ui ��? 1 的部��?

    [umax,in] = max(b(a1));

    [umin,ni] = min(b(a0));

    % a1(in)  得到筛�?�过��? ��? ��?大�?�坐标索��?

    u11(t,:) = [i, umax, a1(in), umin, a0(ni)];  % 里面存放 vi��? 1��? 坐标 对应行的 ��?小�?�和��?大�??

    if length(a1)==1

        s = s + 1;

        loc(s,:) = [i,a1,b(a1)];   % 记录只有��?��? vi=1 ��? 行纵��?

    else
        sss = sss + 1;  % 1的�?�为两个及以��?

        u1{sss,1} = i;  % 记录 行号

        u1{sss,3} = b(a1);  % 第三��? 存储 uik ��? 1 的部分对��? TC 的�??

        u1{sss,2} = a1;  % 单元数组 第二列里面存放中 rik 对应行中 ��?1 的坐��?

    end

    if umin < v1min

      Loc = [i,a0(ni)];  % 记录  ��? ��? ��?

      v1min = umin;  % 更新全局 uik ��?小�??

    end

    t=t+1;
end

% num_row =  cell2mat(u1(:,1));   % 记录行号 形成向量

% num_col =  cell2mat(u1(:,2));   % 记录列�?? 形成向量

% num_value = cell2mat(u1(:,3));  % 记录 相应 TC

% 删除 0 向量
loc(s+1:end,:)=[];
u1(sss+1:end,:)=[];

u0 = zeros(m-v1num,5);  % 此处 重新定义 ��?个类似的单元数组 存放��?0 的部��?

t = v1num + 1;

for i = v0'  % v0 ��?要是行向��?

    [umax,in] = max(TCi(i,:));

    [umin,ni] = min(TCi(i,:));

    u0(t-v1num,:) = [i,umax, in, umin, ni];  % 单元数组 第一列里面存��? vi��? 1的坐��?

    t=t+1;
end

% ��?要给��? vi ��?1的下��?  以及 ��?0 的下��?

% [h l] = max(reshape(TC,m,n), [], 2);   % 给出每一行的��?大�?�列向量h 以及 下标向量l��?  注意1 为每列， 2 为每行�??

model.modelname = 'LB_Inv_UFL';
model.modelsense = 'min';

col = m + m * n;
ncol = col * 2 ;

model.lb  = zeros(ncol, 1);
model.ub  = inf(ncol, 1);

obj = ones(col ,1);

model.obj = [obj; obj];   % norm-1 c-Costs  均为��?

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nplayers, 1)];

% Set data for constraints and matrix

% nrow = m + n + v1num - s + 1; % 前两个约束加起来��? m 个；第三个约束有 n 个；
% 第四个约束有 vi=1 对应行中 rik==1 个数大于 1 的个��? ��? ��?个等��?

nrow = v1num + (m - v1num)*n + n + 1 + (v1num-s) ;

model.A     = sparse(nrow, ncol);

model.sense = [repmat('=', 1, 1); repmat('<', nrow-1, 1)];

% Production constraints   注意限制条件��?要遍��?  这一点非常复��?

model.A(1,:) = x0;

model.rhs(1) = V_UFL-V_0;


%  第一类约��?
for i = 1 : v1num
    % 即需要找��? i ��?在行��? uik ��? 1 的最大�??  % ��?要给出下��?

    inde = zeros(col,1);

    % u1{t,1}()  ��? cell 中取出元��?

    inde((u11(i,1)-1)*n + m + u11(i,3)) = 1;  % max-index

    inde((u11(i,1)-1)*n + m + u11(i,5)) = -1; % min-index

    indexd = [inde; inde*(-1)];

    rhs = (-1)*u11(i,2) + u11(i,4);

    model.A(1 + i, :) = indexd;  % 有问��?

    model.rhs(1 + i) = rhs;

end

%  第二类约��?   % 这里��?  (m-v1num)*s��?
for k= 1:(m - v1num)

  for i=1:s

    inde = zeros(col,1);

    inde(u0(k,1)) = -1; % fi==0的横坐标

    a = v0(k);

    inde((u0(k,1)-1)*n + m + loc(i,2)) = -1;  % 使用 s ��? 纵坐��?

    inde((loc(i,1)-1)*n + m + loc(i,2)) = 1;  % 对应��? rik

    inde(loc(i,1)) = 1;

    indexd = [inde; inde*(-1)];

    rhs = TC((u0(k,1)-1)*n + loc(i,2)) + FC(a) - FC(loc(i,1)) - loc(i,3);

    model.A(v1num + 1 + i + (k-1)*s, :) = indexd;

    model.rhs(v1num + 1 + i + (k-1)*s) = rhs;

  end

end


for k= 1:(m-v1num)  % 这里��?  (m-v1num)*(n-s) ��?

  tt = 1;

  for i = 1 : sss  % 行数

    t=1;

    for j = cell2mat(u1(i,2))

      inde = zeros(col,1);

      inde(u0(k,1)) = -1; % fi==0的横坐标

      a = v0(k);

      inde((u0(k,1)-1)*n + m + j) = -1;  % 使用 s ��? 纵坐��?

      inde((cell2mat(u1(i,1))-1)*n + m + j) = 1;  % 对应��? rik

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

% 第三类约��?   % 这里��?  (s) ��?
for i = 1 : s

   inde = zeros(col,1);

   inde(loc(i,1)) = 1;   % 有问��?

   inde((loc(i,1)-1)*n + m + loc(i,2)) = 1;

   mmi = cell2mat(u1(:,1));

   [mi,mindex] = min(TCi(mmi,loc(i,2))); % ��?小�?�索引即为行��?

   inde((mmi(mindex)-1)*n + m + loc(i,2)) = -1;  %index 有问��?

   indexd = [inde; inde*(-1)];

   rhs = -loc(i,3) - FC(loc(i,1)) + mi;

   model.A(v1num + (m-v1num) * n + 1 + i, :) = indexd;

   model.rhs(v1num + (m-v1num) * n + 1 + i) = rhs;

end

tt = 1;

for i = 1 : sss  % sss行数   % 这里��? (n-s) ��?

  t = 1;

  for j = cell2mat(u1(i,2))

    inde = zeros(col,1);

    inde((cell2mat(u1(i,1))-1)*n + m + j) = 1;   % 对应��? rik

    num_value = cell2mat(u1(i,3));

    mmi = loc(:,1);

    [mi,mindex] = min(TCi(mmi, j));   % ��?个的��?  多个的对应列

    inde((mmi(mindex)-1)*n + m + j) = -1;

    rhs = mi - num_value(t) ;

    indexd = [inde; inde*(-1)];

    model.A(v1num + (m-v1num) * n + 1 + s + tt, :) = indexd;

    model.rhs(v1num + (m-v1num) * n + 1 + s + tt) = rhs;

    t = t + 1;

    tt = tt + 1;

  end

end

% 第四类约��?
for i = 1 : (v1num-s)

  inde = zeros(col,1);

  inde(cell2mat(u1(i,1))) = 1;

  in = cell2mat(u1(i,2))+ m + (cell2mat(u1(i,1))-1)*n; %  位置

  % cell2mat(u1(i,2)) ��? ��? 形成的向��?

  inde(in) = 1;

  new = TCi(loc(1:s,1),cell2mat(u1(i,2)));   % 找到 相应��?

  aa = loc(:,1);

  [mii,ind] = min(sum(new,2));  % ��? ��? 求和  得到 ��?小�?? ��? 行索��?

  inde((aa(ind)-1)*n + m + cell2mat(u1(i,2))) = -1;

  indexd = [inde; inde*(-1)];

  rhs = - sum(cell2mat(u1(i,3))) - FC(cell2mat(u1(i,1))) + mii ;

  model.A(v1num + (m-v1num) * n + 1 + n + i, :) = indexd;

  model.rhs(v1num + (m-v1num) * n + 1 + n + i) = rhs;

end

% gurobi_write(model,'LB_UFL.lp');

result = gurobi(model);

res = result.objval;

end
