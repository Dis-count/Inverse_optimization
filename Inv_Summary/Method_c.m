function Method_c(m,n,k)
  result = zeros(k+1,8);

  for i = 1:k
% k is the number of iteration
    fi = randi([100,200],m,1);
    rij = randi([1,100],m*n,1);

%  rij = randi([1,100],m*n,1);
    ini_sol = feasible_v(m,n);  % Random a feasible solution
%   [opt1,ini_sol] = UFL(fi,rij);  % give the optimal solution
%  s1 s2 s3 are the number of iteration.
    [a,s1] = Cutting(fi,rij,ini_sol);   % The optimal solution

    [b,s2] = Cutting1(fi,rij,ini_sol);  % The UB

    [c,s3] = Cutting2(fi,rij,ini_sol);  % The LB

    result(i,:) = [a,b,c,s1,s2,s3,(b-a)/a,(c-a)/a];

  end

result(k+1,:) = [mean(result(1:k,1)),mean(result(1:k,2)),mean(result(1:k,3)),mean(result(1:k,4)),mean(result(1:k,5)),mean(result(1:k,6)),mean(result(1:k,7)),mean(result(1:k,8))];

filename = ['E:\Files\Matlab\Inv_All\',num2str(m),'by',num2str(n),'.xlsx'];

xlswrite(filename,result);

end
