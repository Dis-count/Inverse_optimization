lambda1 = 4/6;

lambda2 = 20 - 24*lambda1;

f = ones(8,1);

ic = [1:16];

Aeq =[1*lambda1 0 1*lambda2 0 -1*lambda1 0 -1*lambda2 0 ;
    0 1*lambda1 0 1*lambda2 0 -1*lambda1 0 -1*lambda2 ;];

beq =[5 - 6*lambda1 + lambda2 ; 4 - 4*lambda1 - lambda2];

lb = zeros(8,1);

ub = inf(8,1);

[x fval flag] = linprog(f, [], [], Aeq, beq, lb, ub)  % flag 为�??出标�?

% [x, fval] = linprog(f,[],[],A,b,lb,ub)  % 对于 IPU 问题整数规划的解等于线�?�规划的�?
