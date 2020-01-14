function [opt1,opt2] = UFL(fi,rik)

m = length(fi) ;

n = length(rik)/m ;

% Index helper function
%flowidx = @(w, p)  n* p + w;

% Build model
model.modelname = 'I_UFL';
model.modelsense = 'min';

% Set data for variables
ncol = m*n + m ;

model.vtype = 'B';

model.obj   = [fi; rik];

nrow = m*n + n;

model.A     = sparse(nrow, ncol);

model.rhs   = [zeros(m*n, 1); ones(n,1)];

model.sense = [repmat('>', m*n , 1); repmat('=', n , 1)];

for w =1:m
    for p =1:n
        model.A(p+n*(w-1),w) = 1;

        model.A(p+n*(w-1),m+p+(w-1)*n) = -1;
    end
end

for p = 1:n

    model.A(p+m*n, n*(0:(m-1))+m+p) = 1; % �ڶ���Լ��

end

params.outputflag = 0;

res = gurobi(model,params);

opt1 = res.objval;

opt2 = res.x;


end
