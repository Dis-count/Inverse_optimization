function test()
% Copyright 2019, Gurobi Optimization, LLC
% This example formulates and solves the following simple MIP model:
%  maximize
%       5 x +  4 y
%  subject to
%        x + 2 y  <= 7
%        x, y binary

names = {'x'; 'y'};

model.A = sparse([1 2; 2 3]);

model.obj = [5 4];

model.rhs = [7;8];

model.sense = '=';

model.vtype = 'C';

model.modelsense = 'max';

model.varnames = names;

gurobi_write(model, 'Inv.lp');

params.outputflag = 0;

result = gurobi(model, params);

disp(result);

for v=1:length(names)
    fprintf('%s %d\n', names{v}, result.x(v));
end

fprintf('Obj: %e\n', result.objval);
end
