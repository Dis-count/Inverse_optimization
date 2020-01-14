function [a,b,c]= Method_c(m,n,k)
result = zeros(k,6);
  for i = 1:k
% k is the number of iteration
    fi = randi([50,100],m,1);
    rij = randi([10,20],m*n,1);

%  rij = randi([1,100],m*n,1);
    ini_sol = feasible_v(m,n);  % Random a feasible solution
%   [opt1,ini_sol] = UFL(fi,rij);  % give the optimal solution
%  s1 s2 s3 are the number of iteration.
    [a,s1] = Cutting(fi,rij,ini_sol);   % The optimal solution

    [b,s2] = Cutting1(fi,rij,ini_sol);  % The UB

    [c,s3] = Cutting2(fi,rij,ini_sol);  % The LB

    result(k,:) = [a,b,c,s1,s2,s3];

  end

end
