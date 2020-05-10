function  Model2()
%   次梯度方法求解拉格朗日dual
%   max_iter = 100;
%   min_step_size = 0.01;
  x0 = [2;1.5];  % The given solution

  c = [5;4];

  % 松弛约束条件   A x_0 < b
  A = [-1,1;
        6,4;
        1,4;];
% add more constraints

  b = [1;24;9];

  % relax_con = [-1,1,2 ;
  %               3,2,18;];
  iter = 0;
  % con = [1,2,8;]; %  A'\lambda > c

   mu = Min_sub(A,b,c);  %二阶段计算得��? dual

   while iter < 10

     [A0,adjustment] = Model_sub(x0,A,b,c,mu);

     mu = Min_sub(A0,b,c);

     iter = iter +1
%   A0 = reshape(A0,sum_relax,n_row)';
     if abs(b'*mu-c'*x0) < 0.01
       break
     end
   end

  fprintf('A0: %f\n', A0);

  fprintf('mu: %f\n', mu);

  fprintf('iter: %f\n', iter);

  fprintf('adjustment: %f\n', adjustment);

end


function [A0,obj] = Model_sub(x0,A,b,c,mu)
  % min A'-A
  % A'x0 < b ; A'T mu > c ; mu(Ax-b)=0
  % c = [5;4];
  % relax_con = [-1,1,1 ;
  %               6,4,24;];

    model.modelname = 'Model_sub';

    model.modelsense = 'min';

    row0 = length(A(:,1));

    col0 = length(A(1,:));

    ncol = col0*row0;

    nrow = row0 + col0 + row0;

    model.lb = reshape(A',ncol,1);

    model.obj = ones(ncol,1)

    model.vtype = repmat('C',ncol,1);

    AA1 = kron(mu',eye(col0));

    AA = kron(eye(row0),x0');

    B = eye(row0);

    for i = 1:row0

      B(i,:) = B(i,:)*mu(i);

    end

    BB = kron(B,x0');

    model.A = sparse(nrow,ncol);

    for i = 1: row0

      model.A(i,:) = AA(i,:);

    end

    for i = 1: col0

      model.A(row0+i,:) = -AA1(i,:);

    end

    for i = 1: row0

      model.A(row0+col0+i,:) = BB(i,:);

    end

    rhs1 = b;

    rhs2 = -c;

    rhs3 = mu.*b;

    model.rhs = [rhs1;rhs2;rhs3];

    model.sense1 = repmat('<', row0 + col0, 1);

    model.sense2 = repmat('=', row0, 1);

    model.sense =  [model.sense1;model.sense2];

    gurobi_write(model,'model.lp');

    params.outputflag = 0;

    result = gurobi(model, params);

  %  result = gurobi(model)

    opt_x = result.x;  % adjustment

   A0 = reshape(opt_x,row0,col0);

    obj = result.objval - sum(sum(A));

% mu(end)*(mu(1:end-1)'*b-c'*x0)
end


function mu = Min_sub(A,b,c)
  % min bTy
  % c = [5;4];
  %  y > 0
  %  ATy > c
  % Obtain the dual
  A
    model.modelname = 'Min_sub';

    model.modelsense = 'min';

    nrow = length(A(1,:));  %notice that the transposition

    ncol = length(A(:,1));

    model.lb = zeros(ncol,1);

    model.obj = b;

    model.vtype = repmat('C',ncol,1);

    model.A = sparse(nrow,ncol);

    for i = 1: nrow

      model.A(i,:) = A(:,i);

    end

    model.rhs = c;

    model.sense = repmat('>', nrow, 1);

    gurobi_write(model,'min.lp');

    params.outputflag = 0;

    result = gurobi(model, params);

  %  result = gurobi(model)

    mu = result.x;

end
