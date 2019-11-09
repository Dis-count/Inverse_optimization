function [opt_x,cost] = Lang_sub(obj,relax_con,con,mu)

% obj = [5;4];
% relax_con = [6,4,24;
%             -1,1,1 ;];

  model.modelname = 'Sub';

  model.modelsense = 'max';

  nrow = length(con(:,1));

  ncol = length(obj);

  model.lb = zeros(ncol,1);

  model.ub = inf(ncol,1);

for i =1 : length(relax_con(:,1))

  model.obj = obj - relax_con(i,1:ncol)' * mu(i);
%            [5 - 6*mu + nu  ; 4 - 4*mu - nu  ;];

end

  model.vtype = repmat('C',ncol,1);

  model.A = sparse(nrow,ncol);

for i = 1: nrow

  model.A(i,:) = con(i,1:ncol);

end

  model.rhs = con(:,end);

  model.sense = repmat('<', nrow, 1);

  result = gurobi(model);

  opt_x = result.x;

  cost = result.objval + dot(relax_con(:,end),mu) ;

end
