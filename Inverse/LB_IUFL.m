% 这里有四组约��?
% 第一个约束确保行��? r_ik (��?1的最大�??<��?0 的最小�?? )的大小关��?

function res = LB_IUFL(V_UFL,vi,uik,FC,TC)

% V_UFL = 28;   % 给定的最优�?�（目标值）

% vi uik ��? 给定��?优解

% FC TC  ��? 原矩阵�??

m = length(vi);
n = length(uik)/m;

x_0 = [vi; uik];
x0 =[x_0; x_0*(-1)];  %求解��?要的double向量

% 给定的原问题Costs

Costs = [FC; TC];

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

Loc = zeros(1,2);  % 记录 ��?小�?�位��?  约束2用不到了

loc = zeros(v1num,3);  % 记录只有��?个vi��? ��? ��? 坐标 ��? ��?  ��?要删掉剩��? v1num-s ��?
% 约束3用到   约束2也用��?

t = 1;

for i = v1'  % v1 ��?要是行向��? % 循环 每一��?

    a1 = find(ui(i,:) == 1);

    b = TCi(i,:);  % 找出TC中对��? ui ��? 1 的部��?

%     [umax,in] = max(b(a1));
% 
%     a0 = find(ui(i,:) == 0);
% 
%     if length(a0) > 0
% 
%     [umin,ni] = min(b(a0));
% 
%     else
% 
%       umin = 0;
% 
%       ni = 1;
% 
%       a0(ni) = 1;
% 
%     end
%     % a1(in)  得到筛�?�过��? ��? ��?大�?�坐标索��?
% 
%     u11(t,:) = [i, umax, a1(in), umin, a0(ni)];  % 里面存放 vi��? 1��? 坐标 对应行的 ��?小�?�和��?大�??

    if length(a1)==1

        s = s + 1;

        loc(s,:) = [i,a1,b(a1)];   % 记录只有��?��? vi=1 ��? 行纵��?

    else
        sss = sss + 1;  % 1的�?�为两个及以��?

        u1{sss,1} = i;  % 记录 行号

        u1{sss,3} = b(a1);  % 第三��? 存储 uik ��? 1 的部分对��? TC 的�??

        u1{sss,2} = a1;  % 单元数组 第二列里面存放中 rik 对应行中 ��?1 的坐��?

    end

%     if umin < v1min
% 
%       Loc = [i,a0(ni)];  % 记录  ��? ��? ��?
% 
%       v1min = umin;  % 更新全局 uik ��?小�??
% 
%     end

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

if length(v0) >0

for i = v0'  % v0 ��?要是行向��?

    [umax,in] = max(TCi(i,:));

    [umin,ni] = min(TCi(i,:));

    u0(t-v1num,:) = [i,umax, in, umin, ni];  % 单元数组 第一列里面存��? vi��? 1的坐��?

    t=t+1;

end

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

sum1 = (n - 1 + n - s)* s / 2;  % 只有��?��?1 ��?  的�?�约束数

sum2 = 0;

% last_item = length(cell2mat(u1(:,1)));  % 记录 ��?后一��? 总数

if sss > 1

summ = cell2mat(u1(sss,2));

sum3 = length(summ);

for i = (sss-1): -1 : 1

  aa = length(u1(i,2));

  sum2 = sum2 + aa * sum3;   % 记录 总约束数

  sum3 = sum3 + aa;   %  sum3 是一个工具计数变��?

end

end

sum2 = sum2 +sum1;

if s > 0
  nrow = sum2 + (m - v1num)*n + n + 1 + (v1num-s) ;

else

  nrow = sum2 + (m - v1num)*n + n + 1 ;

end

model.A     = sparse(nrow, ncol);

model.sense = [repmat('=', 1, 1); repmat('<', nrow-1, 1)];

% Production constraints   注意限制条件��?要遍��?  这一点非常复��?

model.A(1,:) = x0;

model.rhs(1) = V_UFL-V_0;


%  第一类约��?

[urow, ucol] = find(ui~=0);  % 给出非零向量

urow_col = [urow , ucol];

uu = sortrows(urow_col,2);  % 按列排序

t = 1;  % 记录 约束��?

if (sss ~= 1)||(s ~= 0)

for i = 1 : (n-1)   % ��? 行坐��?

    for j = (i+1) : n

      inde = zeros(col,1);

      % u1(t,1)  cell 中取出元��?
      if uu(i,1) ~= uu(j,1)

      inde((uu(i,1)-1)*n + m + uu(i,2)) = 1;  %  r11

      inde((uu(j,1)-1)*n + m + uu(i,2)) = -1; %  ri1

      inde((uu(j,1)-1)*n + m + uu(j,2)) = 1;  %  ri2

      inde((uu(i,1)-1)*n + m + uu(j,2)) = -1; %  r12

      indexd = [inde; inde*(-1)];

      rhs = - TCi(uu(j,1),uu(j,2)) - TCi(uu(i,1),uu(i,2)) + TCi(uu(j,1),uu(i,2)) + TCi(uu(i,1),uu(j,2));

      model.A(1 + t, :) = indexd;  %

      model.rhs(1 + t) = rhs;

      t = t + 1;

      end

    end
end

end

pt =t;

%  第二类约��?   % 这里��?  (m-v1num)*s��?
if (m - v1num) > 0
  if s > 0

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

      model.A(sum2 + 1 + i + (k-1)*s, :) = indexd;

      model.rhs(sum2 + 1 + i + (k-1)*s) = rhs;

    end

  end

  end

if sss > 0

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


      model.A(sum2 + (m-v1num)*s + 1 + tt + (k-1)*(n-s), :) = indexd;

      model.rhs(sum2 + (m-v1num)*s + 1 + tt + (k-1)*(n-s)) = rhs;

      t = t + 1;

      tt = tt + 1;

    end

  end

end

end

end

% 第三类约��?   % 这里��?  (s) ��?
if s > 0

for i = 1 : s

   inde = zeros(col,1);

   inde(loc(i,1)) = 1;   % 有问��?

   inde((loc(i,1)-1)*n + m + loc(i,2)) = 1;

   mmi = v1;

   mmi(find(mmi==loc(i,1))) = [];

   [mi,mindex] = min(TCi(mmi,loc(i,2))); % ��?小�?�索引即为行��?

   inde((mmi(mindex)-1)*n + m + loc(i,2)) = -1;  %index 有问��?

   indexd = [inde; inde*(-1)];

   rhs = -loc(i,3) - FC(loc(i,1)) + mi;

   model.A(sum2 + (m-v1num) * n + 1 + i, :) = indexd;

   model.rhs(sum2 + (m-v1num) * n + 1 + i) = rhs;

end

end

tt = 1;

if sss > 0

for i = 1 : sss  % sss行数   % 这里��? (n-s) ��?

  t = 1;

  for j = cell2mat(u1(i,2))

    inde = zeros(col,1);

    inde((cell2mat(u1(i,1))-1)*n + m + j) = 1;   % 对应��? rik

    num_value = cell2mat(u1(i,3));

    mmt = v1;

    % mmt(cell2mat(u1(i,1))) = [];  % 注意要删除元�? 而不是位�?

    mmt(find(mmt==cell2mat(u1(i,1)))) = [];

    [mi,mindex] = min(TCi(mmt, j));   % ��?个的��?  多个的对应列   这里是所��?

    inde((mmt(mindex)-1)*n + m + j) = -1;

    rhs1 = mi - num_value(t);

    indexd = [inde; inde*(-1)];

    model.A(sum2 + (m-v1num) * n + 1 + s + tt, :) = indexd;

    model.rhs(sum2 + (m-v1num) * n + 1 + s + tt) = rhs1;

    t = t + 1;

    tt = tt + 1;

  end

end

end

% 第四类约��?
if (v1num-s) > 0

  if s >0

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

    model.A(sum2 + (m-v1num) * n + 1 + n + i, :) = indexd;

    model.rhs(sum2 + (m-v1num) * n + 1 + n + i) = rhs;

  end

  end

end


result1 = gurobi(model);

res = result1.objval;

end
