function opt = Maint(vi0,uij0,c0,I)
% I 表示 restricted set 每一行是给定 feasible solution  行数即可行解个数
% c0 为原设施成本  column vector
% vi0, uij0 表示给定feasible solution [0,1]  Row Vector

m = length(vi0) ;

n = length(uij0)/m ;

model.modelname = 'Main';

model.modelsense = 'min';

ncol = 2*m*n + 2*m ;

model.vtype = 'C';

model.obj   = ones(2*m + 2*m*n,1);

nrow = length(I(:,1));

delta = zeros(nrow,ncol/2);

model.A     = sparse(nrow, ncol);

for p = 1:nrow

    delta(p,:) = I(p,:) - [vi0;uij0]';

    model.A(p,:) = [delta(p,:), -delta(p,:)];

end

model.rhs   = -delta * c0;  % Matrix Multiply

model.sense = repmat('>', nrow , 1);


% Save model
% gurobi_write(model,'Maint.lp');

% Set parameters
% params.method = 2;

% Optimize
% res = gurobi(model, params);
params.outputflag = 0;

res = gurobi(model,params);

opt = res.x;   %  给出 (a,c,b,d)

end
%
% I =[
%      0     1     1     0     0     0     0     0     1     0     1     0     0     1     0     1     0     0     0     0
%      0     0     1     1     0     0     0     0     0     0     0     0     1     0     1     0     0     1     0     1
%      1     1     0     0     1     0     0     1     0     1     1     0     0     0     0     0     0     0     0     0
%      0     1     0     1     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     1
%      0     0     1     0     0     0     0     0     0     0     0     0     1     1     1     1     0     0     0     0
%      1     1     0     0     0     0     0     1     1     1     1     0     0     0     0     0     0     0     0     0
%      0     1     1     0     0     0     0     0     1     0     0     0     0     1     1     1     0     0     0     0
%      0     1     0     1     0     0     0     0     1     1     1     0     0     0     0     0     0     0     0     1];
%
% c0 =[
%     10
%     10
%     10
%     10
%      3
%      3
%     20
%      2
%     20
%      1
%      4
%     20
%      3
%     20
%      3
%      4
%     20
%      2
%     20
%      1]
% vi0 =[0;1;1;0;];
% uij0 = [0; 0; 0; 0; 1; 0; 1; 0; 0; 1; 0; 1; 0; 0; 0; 0;];
