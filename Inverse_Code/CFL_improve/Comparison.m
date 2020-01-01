% This function is used to get excel to show the specific data result.
function Comparison(m,n,k)
% (m,n) 问题规模  k 为调用函数次��? 默认 k=50;

theresult = zeros(k,4);

for j = 1:1

    for i=1:k

      FC = round(rand(1,m)*100)';

      TC = round(rand(1,m*n)*100)';

      [a,b] = fun1(m,n,FC,TC);  % 重复调用

      [c,d] = fun2(m,n,FC,TC);

      theresult(i,:) = [a,b,c,d];

    end

filename = ['E:\Files\Matlab\',num2str(m),'by',num2str(n),'-',num2str(j),'.xlsx'];

% filename = ['E:\Files\Matlab\',num2str(j),'by',num2str(j),'.xlsx'];

xlswrite(filename,theresult);

end

end


function [gap,opt1] = fun1(m,n,FC,TC)

  % ki = round(rand(1,m)*)
  ki = ones(1,m)*100;
  % dj = round(rand(1,n)*)
  dj = ones(1,n);
  % [opt1,opt2] = CFL(FC,TC,ki,dj);  % V_CFL ��? CFL��?优�??

  [opt1,opt2] = UFL(FC,TC);

  V_UFL = opt1;

  vi = opt2(1:m);

  uik = opt2(m+1:end);
  % Build model
  model.modelname = 'H_UFL';
  model.modelsense = 'min';

  % Set data for variables
  ncol = 5*m*n + 4*m + n;

  % 先试试变量大于零的情��?
  model.lb    = zeros(ncol, 1);
  model.ub    = [inf(ncol, 1)];
  model.obj   = [zeros(n+m+3*m*n,1); ones(2*m + 2*m*n,1); zeros(m,1);];
  % model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nWarehouses, 1)];
  %
  % for p = 1:nPlants
  %     model.varnames{p} = sprintf('Open%d', p);
  % end
  %
  % for w = 1:nWarehouses
  %     for p = 1:nPlants
  %         v = flowidx(w, p);
  %         model.varnames{v} = sprintf('Trans%d,%d', w, p);
  %     end
  % end

  % Set data for constraints and matrix
  nrow = 2*m*n+2*m+2;

  model.A     = sparse(nrow, ncol);

  model.rhs   = [V_UFL; zeros(m + m*n, 1); V_UFL; FC; TC];
  model.sense = [repmat('>', 1, 1); repmat('=', 2*m*n + 2*m + 1, 1)];

  model.A(1,1:n) = 1;  % 第一类约��?

  for p = 1:m
      for w = 1:n
          model.A(p+1, n*p+w) = 1;
      end
      model.A(p+1, n+n*m+p) = -1;
      model.A(p+1, 5*m*n + 3*m +n+p) = ki(p)
  %    model.constrnames{p} = sprintf('Capacity%d', p);
  end

  % 第二类约��?

  for p = 1:m
      for w = 1:n
          model.A((p-1)*n+w+m+1,[w,m*n+m+p*n+w]) = 1;

          model.A((p-1)*n+w+m+1,[p*n+w,2*m*n+m+p*n+w]) = -1;
          model.A((p-1)*n+w+m+1, 5*m*n + 3*m +n+p) = -dj(w)
      end
  end   % 第三类约��?

  for p = 1:m
      for w = 1:n
      model.A(m*n+m+2,m+2*m*n+p*n+w) = uik((p-1)*n+w);

      end
      model.A(m*n+m+2,n+m*n+p) = vi(p);
  %    model.constrnames{nPlants+w} = sprintf('Demand%d', w);
  end   % 第四个约��?

  for p = 1:m
      model.A(m*n+m+2+p,n+m+3*m*n+p) = -1;
      model.A(m*n+m+2+p,[n+m*n+p,n+2*m+3*m*n+p]) = 1;   % 保持右侧约束为正的fi 下同
  end   % 第五个约��?


  for p = 1:m
      for w = 1:n
          model.A((p-1)*n+m*n+2*m+2+w, 3*m+3*m*n+p*n+w) = -1;

          model.A((p-1)*n+m*n+2*m+2+w, [m+2*m*n+p*n+w,3*m+4*m*n+p*n+w]) = 1;
      end
  end  % 第六个约��?


  % Save model
  % gurobi_write(model,'UFL.lp');
  % Guess at the starting point: close the plant with the highest fixed
  % costs; open all others first open all plants

  % model.start = [ones(m, 1); inf(nPlants * nWarehouses, 1)];
  % [~, idx] = max(FixedCosts);
  % model.start(idx) = 0;

  % Set parameters
  % params.method = 2;

  % Optimize
  % res = gurobi(model, params);
  params.outputflag = 0;

  result = gurobi(model,params);

  gap = result.objval;

  end


  function [gap,opt1] = fun2(m,n,FC,TC)  % ����ԭUFL������ֵ �� �� LP �����õ��� gap
  % ���� Facility ���� m  Player ���� n

  % mul Ϊ ��ʩ�ɱ�����

  % ��Ӧ������������ Ϊ ��2mn+2m+2����  �����У�5mn+3m+n) ��

  % rik = randint(1,m,[1000 2000])';
  % rik = fun3(m,n,round(m*n*0.25));


  [opt1,opt2] = UFL(FC,TC);  % V_UFL Ϊ UFL����ֵ
  V_UFL = opt1;

  % fi  = [10; 10; 10; 10];
  % �������� facility cost

  % M = 30;    % define M= a bigger integer.


  % �������� transportation cost
  % rik  = [ 3; 3; M; 2;
  %          M; 1; 4; M;
  %          3; M; 3; 4;
  %          M; 2; M; 1;];


  vi = opt2(1:m);

  uik = opt2(m+1:end);

  % ����һ�� ����Ϊ ��5mn+3m+n) ��������ÿһ���������� ���� һ��������

  % Build model
  model.modelname = 'UFL';
  model.modelsense = 'min';

  % Set data for variables
  ncol = 5*m*n + 3*m + n;

  % �����Ա���������������
  model.lb    = zeros(ncol, 1);
  model.ub    = inf(ncol, 1);
  model.obj   = [zeros(n+m+3*m*n,1); ones(2*m + 2*m*n,1); ];

  % Set data for constraints and matrix
  nrow = 2*m*n+2*m+2;

  model.A     = sparse(nrow, ncol);


  model.rhs   = [V_UFL; zeros(m + m*n, 1); V_UFL; FC; TC];
  model.sense = [repmat('>', 1, 1); repmat('=', 2*m*n + 2*m + 1, 1)];


  model.A(1,1:n) = 1;  % ��һ��Լ��

  for p = 1:m
      for w = 1:n
          model.A(p+1, n*p+w) = 1;
      end
      model.A(p+1, n+n*m+p) = -1;
  %    model.constrnames{p} = sprintf('Capacity%d', p);
  end

  % �ڶ���Լ��

  for p = 1:m
      for w = 1:n
          model.A((p-1)*n+w+m+1,[w,m*n+m+p*n+w]) = 1;

          model.A((p-1)*n+w+m+1,[p*n+w,2*m*n+m+p*n+w]) = -1;
      end
  end   % ������Լ��

  for p = 1:m
      for w = 1:n
      model.A(m*n+m+2,m+2*m*n+p*n+w) = uik((p-1)*n+w);

      end
      model.A(m*n+m+2,n+m*n+p) = vi(p);
  %    model.constrnames{nPlants+w} = sprintf('Demand%d', w);
  end   % ���ĸ�Լ��


  for p = 1:m
      model.A(m*n+m+2+p,n+m+3*m*n+p) = -1;
      model.A(m*n+m+2+p,[n+m*n+p,n+2*m+3*m*n+p]) = 1;   % �����Ҳ�Լ��Ϊ����fi ��ͬ
  end   % ������Լ��


  for p = 1:m
      for w = 1:n
          model.A((p-1)*n+m*n+2*m+2+w, 3*m+3*m*n+p*n+w) = -1;

          model.A((p-1)*n+m*n+2*m+2+w, [m+2*m*n+p*n+w,3*m+4*m*n+p*n+w]) = 1;
      end
  end  % ������Լ��

  % Save model
  % gurobi_write(model,'UFL.lp');

  % Optimize
  res = gurobi(model);

  gap = res.objval;

  end
