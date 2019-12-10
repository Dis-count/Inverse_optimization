function [opt_x,cost] = Sub_I(obj,relax_con,con,mu)

% obj = [16;10;0;4];
% relax_con = [8;2;1;4];

  model.modelname = 'Sub';

  model.modelsense = 'max';

  nrow = length(con(:,1));

  ncol = length(obj);

  model.lb = zeros(ncol,1);

  model.ub = inf(ncol,1);

  model.obj = obj - relax_con *mu;
%            [16 - 8*mu; 10 - 2*mu; -mu; 4 - 4*mu];

  model.vtype = repmat('B',ncol,1);

  model.A = sparse(nrow,ncol);

  model.A(1,:) = con(1,:);

  model.A(2,:) = con(2,:);

  model.rhs = [1,1];

  model.sense = repmat('<', nrow, 1);

  result = gurobi(model);

  opt_x = result.x;

  cost = result.objval + 10*mu ;

end
