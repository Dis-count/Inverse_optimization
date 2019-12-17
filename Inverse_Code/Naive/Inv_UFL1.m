function Inv_UFL1()    % 有四重循环，且只针对m=n=4的情况，需要进一步优化
%  x_0 为给定服务对象以及是否开启 Facility  v_0 为目标值
fixed = 30;

v_i = [0;0;1;1];

u_ij = [0;0;0;0;
        0;0;0;0;
        1;0;1;0;
        0;1;0;1;];

x_0 = [v_i;u_ij];

x0 =[x_0;x_0*(-1)];   % double the matrix

% define primitive data
m = 4;   % Facility
n = 4;   % Players

% Warehouse demand in thousands of units
% Demand      = [15; 18; 14; 20];

% Plant capacity in thousands of units
% Capacity    = [20; 22; 17; 19; 18];
% Fixed costs for each plant

FC  = [10; 10; 10; 10];
% Transportation costs per thousand units
M = 30;
TC  = [ 3; 3; M; 2;
         M; 1; 4; M;
         3; M; 3; 4;
         M; 2; M; 1;];
Costs =[FC;TC];    % Cost 为给定

V_0 = x_0'*Costs;    %  原最优解


% Index helper function
% flowidx = @(w, p) nPlants * w + p;

% Build model
model.modelname = 'Inv_facility';
model.modelsense = 'min';

% Set data for variables
ncol = (m + m * n)*2 ;
model.lb    = zeros(ncol, 1);
model.ub    = [inf(ncol, 1)];
obj = [ones(6,1);0;1;0;1;1;0;1;0;1;1;0;1;0;1];

model.obj = [obj; obj];   % norm-1 c-Costs  均为正
% model.vtype = [repmat('B', nPlants, 1); repmat('C', nPlants * nplayers, 1)];

% for p = 1:nPlants
%     model.varnames{p} = sprintf('Open%d', p);
% end
%
% for w = 1:nplayers
%     for p = 1:nPlants
%         v = flowidx(w, p);
%         model.varnames{v} = sprintf('Trans%d,%d', w, p);
%     end
% end

% Set data for constraints and matrix
nrow = m^n + 1;
model.A     = sparse(nrow, ncol);

model.sense = [repmat('=', 1, 1); repmat('>', m^n, 1)];

% Production constraints   注意限制条件需要遍历  这一点非常复杂
% for p = 1:nPlants
%     for w = 1:nplayers
%         model.A(p, p) = -Capacity(p);
%         model.A(p, flowidx(w, p)) = 1.0;
%     end
%     model.constrnames{p} = sprintf('Capacity%d', p);
% end

% Demand constraints
% for w = 1:nplayers
%     for p = 1:nPlants
%         model.A(nPlants+w, flowidx(w, p)) = 1.0;
%     end
%     model.constrnames{nPlants+w} = sprintf('Demand%d', w);
% end

model.A(1,:) = x0;
model.rhs(1) = fixed-29;
a = eye(m);
% f = zeros(m^n,m+m*n);
s = 2;
for i = 1:m
    for j = 1:m
        for k = 1:m
            for h = 1:m

                b = [a(:,i);a(:,j);a(:,k);a(:,h)];
                c = reshape(b,m,m);
                f = reshape(c',m*m,1);

                d = [(i==1)|(j==1)|(k==1)|(h==1);
                     (i==2)|(j==2)|(k==2)|(h==2);
                     (i==3)|(j==3)|(k==3)|(h==3);
                     (i==4)|(j==4)|(k==4)|(h==4);];
                e = [d;f];   % 转为行向量
                ee = [e;e*(-1)];
                model.A(s,:) = ee;
                model.rhs(s) = fixed - e'* Costs;
                % f(s,:) = e;
                s = s+1;

            end
        end
    end
end
%  循环结束后生成限制条件的所有向量


% Save model
% gurobi_write(model,'Inv_UFL.lp');

% Guess at the starting point: close the plant with the highest fixed
% costs; open all others first open all plants
% model.start = [ones(nPlants, 1); inf(nPlants * nplayers, 1)];
% [~, idx] = max(FC);
% model.start(idx) = 0;



% Optimize
res = gurobi(model);

% Print solution
% if strcmp(res.status, 'OPTIMAL')
%     fprintf('\nTotal Costs: %g\n', res.objval);
%     fprintf('solution:\n');
%     for p = 1:nPlants
%         if res.x(p) > 0.99
%             fprintf('Plant %d open:\n', p);
%         end
%         for w = 1:nplayers
%             if res.x(flowidx(w, p)) > 0.0001
%                 fprintf('  Transport %g units to warehouse %d\n', res.x(flowidx(w, p)), w);
%             end
%         end
%     end
% else
%     fprintf('\n No solution\n');
% end

end
