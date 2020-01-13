function x = Cutting2(vi0,uij0,fi,rij)
% This one is used to obtain the lower bound by local search.
% vi0 = [0;1;1;0;];
% uij0 = [0; 0; 0; 0; 1; 0; 1; 0; 0; 1; 0; 1; 0; 0; 0; 0;];
% fi = [10; 10; 10; 10;];
% vi0, uij0 表示给定 feasible solution [0,1]  row vector
% x 表示主问题给出的 (a,c,b,d)  2mn+2m
% c0 = (fi;rij) 为原设施成本  column vector
%  You must notice that the row and column vectors!!!!
m = length(fi);
n = length(rij)/m;

s = 0; % 计数
mycost = -0.5;
% 初始化限制集
I = [vi0;uij0]';

c0 = [fi;rij];

while mycost < -0.00001
  x = Main(vi0,uij0,c0,I);

  c = (x(1:m*n+m)-x(m*n+m+1:2*m*n+2*m)+[fi;rij]);

  opt2 = dot(c,[vi0;uij0]);

  opt1 = [vi0;uij0];

  [mycost,opt_sol] = LocalS(opt1,opt2,c,m,n); % Local search give the neighborhood solution

  vi0 = opt_sol(1:m);
  uij0 = opt_sol(m+1:end);

  I = [I;opt_sol'];

  s = s + 1;

  if s > 80
    break
  end

end
s
% function end: 'myFunction'
end


function opt = Main(vi0,uij0,c0,I)
% I 表示 restricted set 每一行是给定 feasible solution  行数即可行解个数
% c0 为原设施成本  column vector
% vi0, uij0 表示给定feasible solution [0,1]  Row Vector

m = length(vi0) ;

n = length(uij0)/m ;

model.modelname = 'Main';

model.modelsense = 'min';

ncol = 2*m*n + 2*m ;

model.vtype = 'C';

model.obj   = ones(2*m + 2*m*n,1);

nrow = length(I(:,1));

delta = zeros(nrow,ncol/2);

model.A     = sparse(nrow, ncol);

for p = 1:nrow

    delta(p,:) = I(p,:) - [vi0;uij0]';

    model.A(p,:) = [delta(p,:), -delta(p,:)];

end

model.rhs   = -delta * c0;  % Matrix Multiply

model.sense = repmat('>', nrow , 1);


% Save model
% gurobi_write(model,'Main.lp');

% Set parameters
% params.method = 2;

% Optimize
% res = gurobi(model, params);
params.outputflag = 0;

res = gurobi(model,params);

opt = res.x;   %  给出 (a,c,b,d)

end


function [mycost,opt_sol] = LocalS(opt1,opt2,c,m,n)
% opt1 is a suboptimal solution
% c0 = (fi;rij) 为原设施成本  Column vector
% opt is a local optimal solution
% c is the objective parameters

fi = c(1:m);

rij = c(m+1:end);

opt0 = dot(c,opt1) ;  % Initial feasible cost

% Define the number of opening and not opening facilities
facility = opt1(1:m);  % Get the first m vectors.
trans = reshape(rij,n,m)'; % Reshape to a m*n cost matrix

open1 = sum(facility);
ind_open = find(facility==1); % Record the initial index of the corresponding facilities.

unopen = m - open1;
ind_unopen = setdiff(1:m,ind_open)'; % Obtain the Set Difference

optcost = opt0;

opt_sol = opt1;

for i = 1:2*(m+n)  % Control the number of iteration.

  if unopen > 0
    % Define the add-move   p kinds of local search
    [opt_sol1,cost1] = add0(fi,trans,ind_open,ind_unopen);
  end
  % Define the swap-move   p(m-p) kinds of local search
  if unopen > 0

    [opt_sol2,cost2] = swap(fi,trans,ind_open,ind_unopen);
  end

% Define the remove   m-p  kinds of local search
  if open1 > 1

    [opt_sol3,cost3] = remove0(fi,trans,ind_open);
  end

  [opt_cost,opt_ind] = min([cost1,cost2,cost3]);

  optsol = [opt_sol1,opt_sol2,opt_sol3];

  if opt_cost < optcost

    optcost = opt_cost; % obtain the suboptimal cost

    opt_sol = optsol(:,opt_ind);

    opt_fac = opt_sol(1:m);

    ind_open = find(opt_fac==1); % The index of opening

    open1  = sum(opt_fac);

    unopen = m - open1;

    ind_unopen = setdiff(1:m,ind_open)';

  end

end

mycost = optcost - opt2;

end
