function [mycost,opt_sol] = Local(opt1,opt2,c,m,n)
% opt1 is a suboptimal solution
% c0 = (fi;rij) 为原设施成本  Column vector
% opt is a local optimal solution
% c is the objective parameters

fi = c(1:m);

rij = c(m+1:end);

opt0 = dot(c,opt1) ;  % Initial feasible cost

% Define the number of opening and not opening facilities
facility = opt1(1:m);  % Get the first m vectors.
trans = reshape(rij,n,m)'; % Reshape to a m*n cost matrix

open1 = sum(facility);
ind_open = find(facility==1); % Record the initial index of the corresponding facilities.

unopen = m - open1;
ind_unopen = setdiff(1:m,ind_open)'; % Obtain the Set Difference

optcost = opt0;

for i = 1:2*(m+n)  % Control the number of iteration.

  if unopen > 0
    % Define the add-move   p kinds of local search
    [opt_sol1,cost1] = add0(fi,trans,ind_open,ind_unopen);
  end
  % Define the swap-move   p(m-p) kinds of local search
  if unopen > 0

    [opt_sol2,cost2] = swap(fi,trans,ind_open,ind_unopen);
  end

% Define the remove   m-p  kinds of local search
  if open1 > 1

    [opt_sol3,cost3] = remove0(fi,trans,ind_open);
  end

  [opt_cost,opt_ind] = min([cost1,cost2,cost3]);

  optsol = [opt_sol1,opt_sol2,opt_sol3];

  if opt_cost < optcost

    optcost = opt_cost; % obtain the suboptimal cost

    opt_sol = optsol(:,opt_ind);

    opt_fac = opt_sol(1:m);

    ind_open = find(opt_fac==1); % The index of opening

    open1  = sum(opt_fac);

    unopen = m - open1;

    ind_unopen = setdiff(1:m,ind_open)';

  end

end

mycost = optcost - opt2;

end


% opt1=[
%       0
%       1
%       1
%       0
%       0
%       0
%       0
%       0
%       1
%       0
%       1
%       0
%       1
%       1
%       0
%       0
%       0
%       0
%       0
%       0];

% opt2 =0;
%
% M =20;
%
% c =[10; 10; 10; 10;
%     3; 3; M; 2;
%     M; 1; 4; M;
%     3; M; 3; 4;
%     M; 2; M; 1;];
%  m=4;
%  n=4;
%
