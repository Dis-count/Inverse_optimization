% This function is used to get excel to show the specific data result.
function Ran_CFL(m,n,k)
% (m,n) 问题规模  k 为调用函数次��? 默认 k=50;

theresult=zeros(k,2);

for j = 2:10

    for i=1:k

        [a,b] = Fun1(m,n,j);  % 重复调用
        theresult(i,:) = [a,b];

    end

filename = ['E:\Files\Matlab\',num2str(m),'by',num2str(n),'-',num2str(j),'.xlsx'];

% filename = ['E:\Files\Matlab\',num2str(j),'by',num2str(j),'.xlsx'];

xlswrite(filename,theresult);

end

end

function [gap,opt1] = Fun1(m,n,mul)

  FC = round(rand(1,m)*10*mul)';
  % 可以调整 fi ��? rik 之间的比��? 查看gap变化
  TC = round(rand(1,m*n)*10)';

  % ki = round(rand(1,m)*)
  ki = ones(1,m)*100;
  % dj = round(rand(1,n)*)
  dj = ones(1,n);

  [opt1,opt2] = CFL(FC,TC,ki,dj);  % V_CFL ��? CFL��?优�??

  V_UFL = opt1;

  vi = opt2(1:m);

  uik = opt2(m+1:end);

  % Build model
  model.modelname = 'CFL';
  model.modelsense = 'min';

  % Set the number of variables

  ncol = 4*m*n + 4*m + n;

  % variable > 0

  model.lb    = zeros(ncol, 1);
  model.ub    = [inf(ncol, 1)];
  model.obj   = [zeros(n+2*m+2*m*n,1); ones(2*m + 2*m*n,1);];

  % model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nWarehouses, 1)];

  % Set the number of constraints.

  nrow = 2*m*n+2*m+2;

  model.A     = sparse(nrow, ncol);

  model.rhs   = [V_UFL; zeros(m + m*n, 1); V_UFL; FC; TC];

  model.sense = [repmat('>', 1, 1); repmat('=', 2*m*n + 2*m + 1, 1)];

  model.A(1,1:n) = 1;  % The first class of constraints

  for p = 1:m

      model.A(p+1, n+p) = ki(p);

      model.A(p+1, n+m+p) = -1;

  end

  % The Second class of constraints

  for p = 1:m
      for w = 1:n
          model.A((p-1)*n+w+m+1,[w,2*m+p*n+w]) = 1;

          model.A((p-1)*n+w+m+1,m*n+2*m+p*n+w) = -1;

          model.A((p-1)*n+w+m+1,n+p) = -dj(w);

      end
  end   % The Third class of constraints


  for p = 1:m
      for w = 1:n

      model.A(m*n+m+2,2*m+m*n+p*n+w) = uik((p-1)*n+w);

      end

      model.A(m*n+m+2,n+m+p) = vi(p);

  end   % The Fourth class of constraints


  for p = 1:m
      model.A(m*n+m+2+p,n+2*m+2*m*n+p) = -1;

      model.A(m*n+m+2+p,[n+m+p,n+3*m+2*m*n+p]) = 1;   % 保持右侧约束为正的fi 下同
  end   % The Fifth class of constraints


  for p = 1:m
      for w = 1:n
          model.A((p-1)*n+m*n+2*m+2+w, 4*m+2*m*n+p*n+w) = -1;

          model.A((p-1)*n+m*n+2*m+2+w, [2*m+m*n+p*n+w,4*m+3*m*n+p*n+w]) = 1;
      end
  end  % The Sixth class of constraints

  % Save model
  % gurobi_write(model,'UFL.lp');

  % Optimize
  % res = gurobi(model, params);

  params.outputflag = 0;

  result = gurobi(model,params);

  gap = result.objval;

% function end: 'myFunction'
end
