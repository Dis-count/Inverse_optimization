function x_0 = Sub_Max(A,b,c)

% A = [6 4.8; -1 1.2];
% b =[24;1];
% c =[5;4];

model.modelname = 'Inv_Master_Sub';

model.modelsense = 'max';

[ncol, nrow] = size(A);  % length(y) = length(c) = ncol;

% Set data for variables

model.lb    = zeros(ncol, 1);

model.ub    = [inf(ncol, 1)];

model.obj = b;

% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants* nplayers, 1)];

% Set data for constraints and matrix

%  model.A     = sparse(nrow, ncol);

model.sense = repmat('=', ncol, 1);

model.A = sparse(A');

model.rhs = c;

gurobi_write(model, 'Inv_Master_Sub.lp');

result = gurobi(model);

% res = result.objval;

x_0 = result.x;

% function end: 'myFunction'

end
