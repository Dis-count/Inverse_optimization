best_ub = 1e10;

best_lb = 0;

best_mu = 0;

best_sl = zeros(1,4);



function  = subsolve(min_step_size, max_iter)
%   次梯度方法求解拉格朗日对偶
   iter = 0;

   non_improve = 0;

   max_non_improve = 3;

   lambda = 2;

   step_size = 1;

   mu = 0;   % 初始化拉格朗日乘子

   sp(mu)   % 松弛第一个约束条件

while iter < max_iter

  sp.chang(mu);

  if  sp.solve() == false

    disp('solve wrong!');

  end

% 更新上界
  if sp.opt_cost < best_ub

    best_ub = sp.opt_cost;

    best_mu = mu;

    for i = 0: length(best_sl)-1

      best_sl(i) = sp.opt_x(i);

    end

    non_improve = 0;

  else

    non_improve = non_improve + 1;

  end

  fprintf('iter: %d\n', iter);

  fprintf('best lb: %f ', best_lb);

  fprintf('best ub: %f ', best_ub);

  fprintf('mu: %f ', mu);

  subgradient = 8*sp.opt_x[0] + 2*sp.opt_x[1] + sp.opt_x[2] + 4*sp.opt_x[3] - 10;

  mu = max(0, mu + step_size * subgradient);

% 满足原问题约束的可行解可以作为原问题的下界

  if subgradient <= 0

    current_lb = 16*sp.opt_x[0] + 10*sp.opt_x[1] + 4*sp.opt_x[3];

    if current_lb > best_lb

      best_lb = current_lb;

    end

  end

% 上界未更新达到一定次数

  if non_improve >= max_non_improve

    lambda = lambda/2;

    non_improve = 0;

  end

  mydist = subgradient^2;

  % 迭代停止条件2和3

  if (dist <= 0)||(best_lb >= best_ub-0.0000001)

    break;

  end

  step_size = lambda * (sp.opt_cost - best_lb)/mydist;

  % 迭代停止条件4

  if step_size < min_step_size

    break;

  end

end
