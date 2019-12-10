function [opt_x,cost] = Sub_L(obj,relax_con,con,mu)

% obj = [5;4];
% Notice here that relax_con have several constraints

% relax_con = [6,4,24;];

  model.modelname = 'Sub';

  model.modelsense = 'max';

  nrow = length(con(:,1));

  ncol = length(obj);

  model.lb = zeros(ncol,1);

  model.ub = inf(ncol,1);

  model.obj = obj - relax_con(1:ncol)' *mu;
%            [16 - 8*mu; 10 - 2*mu; -mu; 4 - 4*mu];

  model.vtype = repmat('C',ncol,1);

  model.A = sparse(nrow,ncol);

for i = 1: nrow

  model.A(i,:) = con(i,1:ncol);

end

  model.rhs = con(:,end);

  model.sense = repmat('<', nrow, 1);

  result = gurobi(model);

  opt_x = result.x;

  cost = result.objval + relax_con(end)*mu ;

end
