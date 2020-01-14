function sol = feasible_v(m,n)

  vi0 = zeros(m,1);
  uij0 = zeros(m*n,1); % The initial values are all zeros.
% This function is used to obtain the feasible solution.
% m is the number of facility
% n is the number of custome  n>=m
  k = randi(m);   % Select the number of opening facility

  open_faci = randperm(m,k);  % Select k numbers from [1,n] unrepeatedly

  vi0(open_faci) = 1;

  open_fac = sort(open_faci); % facilitate to obtain the sequence

  a = randi(k,1,n-k); % Return the repeated number [1,k] select n-k times

  b = [a,1:k];  % The selected k facilities must have at least one custome.

  % In fact, it is a process of grouping.

  b = b(randperm(length(b)));  % 随机打乱一个数组的顺序

  open_seq = open_fac(b);  % Obtain the real row sequence
  % The corresponding column sequence is 1:n

  trans_index = [(open_seq-1)*n]+[1:n];   % The whole sequence can be obtained
  uij0(trans_index) = 1;

  sol = [vi0;uij0];

end
