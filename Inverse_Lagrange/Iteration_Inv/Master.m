%   The model   min |A'-A|   A 是原矩阵
%          s.t. (A')y   =   c
%               c x^0 \geq  b y

%  c = [5;4];
%  A = [6 4;-1 1];
%  b = [24;1];
%  y = [\lambda1;\lambda2];  optimal = [5/6;0];
%  We can set anyone at first
% x0 = [2;2.5];

function x_0 = Master(A,b,c,y,x0)

% x0 is feasible \to optimal
% b is the rhs.
% c is the objective coefficient.
% A is the System Matrix.
% y is the dual decision variable. which satisfy the cx^0 = by
% row vector

model.modelname = 'Inv_Master';

model.modelsense = 'min';

[ncol, nrow] = size(A);  % 注意约束��? 是转��?

length_A = ncol * nrow;
% Set data for variables

model.lb  = zeros(length_A*2, 1);

model.ub  = [inf(length_A*2, 1)];

model.obj = ones(length_A*2, 1);   % norm-1 c-Costs  均为��?

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants* nplayers, 1)];

% Set data for constraints and matrix

model.A     = sparse(nrow, length_A*2);

model.sense = repmat('=', nrow, 1);

for i = 1: nrow

    model.A(i,1:length_A) = [zeros(1,(i-1)*nrow),y',zeros(1,length_A-nrow*i)];

    model.A(i,length_A+1:2*length_A) = -[zeros(1,(i-1)*nrow),y',zeros(1,length_A-nrow*i)];
end

model.rhs = c - A'*y;

gurobi_write(model, 'Inv_Master.lp');

% Guess at the starting point: close the plant with the highest fixed
% costs; open all others first open all plants
% model.start = [ones(nPlants, 1); inf(nPlants * nplayers, 1)];
% [~, idx] = max(FC);
% model.start(idx) = 0;

% Optimize

result = gurobi(model);

res = result.objval;

x_0 = result.x;

% function end: 'myFunction'

end
