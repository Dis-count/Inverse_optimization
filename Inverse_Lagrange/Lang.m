function  Lang(min_step_size, max_iter)
%   次梯度方法求解拉格朗日对偶
  best_ub = 1e6;

  best_lb = -1;

%   max_iter = 100;
%   min_step_size = 0.01;

  obj = [5;4];

  iter = 1;

  non_improve = 0;

  max_non_improve = 3;

  % 松弛两个约束条件
  % relax_con = [-1,1,1;
  %               6,4,24;];

  relax_con = [-1,1,2 ;
                3,2,18;];

  con = [1,2,8;];

  sum_relax = length(relax_con(:,1));

%  lambda = ones(sum_relax,1);
  lambda = 0.5;   % 步长

  subgradient = zeros(sum_relax,1);

  step_size = 1;

  mu = zeros(sum_relax,1);   % 初始化拉格朗日乘子

while iter < max_iter

  [opt_x,opt_cost] = Lang_sub(obj,relax_con,con,mu);

  % if  sp.solve() == false
  %
  %   disp('solve wrong!');
  %
  % end

opt_x

% 更新上界
  if opt_cost < best_ub

    best_ub = opt_cost;

    non_improve = 0;

  else

    non_improve = non_improve + 1;

  end

% Notice here that the subgradient is a vector.
  for i = 1 : sum_relax

    subgradient(i) = dot(opt_x,relax_con(i,1:(end-1))) - relax_con(i,end);

    mu(i) = max(0, (mu(i) + step_size * subgradient(i)));

  end

% 满足原问题约束的可行解可以作为原问题的下��?

  if all(subgradient <= 0)   % 如果 subgradient <0 说明满足原问题约��?

    current_lb = dot(opt_x,obj);

    if current_lb > best_lb

      best_lb = current_lb;

    end

  end

  fprintf('iter: %d\n', iter);

  fprintf('best lb: %f ', best_lb);

  fprintf('best ub: %f\n ', best_ub);

  mu
%  fprintf('mu: %f %f\n\n', mu);

% 上界未更新达到一定次��?

  if non_improve >= max_non_improve

    lambda = lambda/2;

    non_improve = 0;

  end

  mydist = norm(subgradient);    % Obtain the norm of subgradient

  % 迭代停止条件2��?3

  if (any(mydist) <= 0)||(best_lb >= best_ub-0.00001)

    break;

  end

  step_size = lambda * (opt_cost - best_lb)/mydist;

  % 迭代停止条件4

  if any(step_size) < min_step_size

    break;

  end

iter = iter + 1;

end
