
% The original simple inverse is stated as follows
% min |a_1 -1| + |b_1-1| + |a_2 -6| + |b_2 -4|
% -2a_1 + 2.5b_1 = 1
% 2 a_2 + 2.5b_2 = 24

% which can be transfer to the following formation according to the L-1 norm
% double the variable

% This function is used to calculate the solution under the L-1 norm.

% A = [2;2.5;2;2.5];
% c_0 = [-1;1;6;4];
% In fact the corresponding coefficient is A = [-2;2.5;0;0; 2;2.5;0;0]
% b = [1;24];

function x_0 = Simple_test(A,b,c_0)
% Notice that A is the original left-side matrix and it's ordered by the row.
% b is the corresponding right-side matrix
% c_0 is the objective coefficient don't need the zero vectors.
% so the basic model is
%      min |x-c_0|
% s.t.  Ax \geq b
% So it's obvious that b need to be changed

% double the matrix

%  给定的最优解向量 对应的原成本矩阵的花费�??

% Build model
model.modelname = 'Inv_Linear_Production';

model.modelsense = 'min';

nrow = length(b);

Acol = length(A);

col = length(A)/nrow;

% Set data for variables
ncol = Acol*2 ;

A_0 = reshape(A,Acol/nrow,nrow);
A_0 = A_0';

cc = reshape(c_0,Acol/nrow,nrow);
cc = cc';

model.lb    = zeros(ncol, 1);

model.ub    = [inf(ncol, 1)];

obj = ones(length(c_0),1);

model.obj = [obj; obj];   % norm-1 c-Costs  均为��?

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants* nplayers, 1)];

% Set data for constraints and matrix

model.A     = sparse(nrow, ncol);

model.sense = repmat('=', nrow, 1);

% Production constraints

for i = 1: nrow

    model.A(i,:) = [zeros((i-1)*col,1); A_0(i,:)';zeros((nrow-i)*col,1); zeros((i-1)*col,1); A_0(i,:)'*(-1); zeros((nrow-i)*col,1)];

    model.rhs(i) = b(i) - dot(cc(i,:), A_0(i,:)) ;

end

gurobi_write(model,'Inv_Linear.lp');

% Guess at the starting point: close the plant with the highest fixed
% costs; open all others first open all plants
% model.start = [ones(nPlants, 1); inf(nPlants * nplayers, 1)];
% [~, idx] = max(FC);
% model.start(idx) = 0;

% Optimize

result = gurobi(model);

res = result.objval;
x_0 = result.x;


end
