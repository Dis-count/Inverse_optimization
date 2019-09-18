% 这里有四组约束
% 第一个约束确保行内 r_ik 的大小关系


V_UFL = 28;   % 给定的最优值（目标值）

vi = [1 ; 1 ; 0;];

uik = [ 1; 0 ; 1 ; 
        0; 1 ; 0 ;  
        0; 0 ; 0 ;];

m = length(vi);  
n = length(uik)/m;    
    
x_0 = [vi;uik];
x0 =[x_0;x_0*(-1)];  %求解需要的double向量


% 给定的原问题Costs

FC    = [5; 6; 7;];

TC    = [11; 4; 8; 
         5; 7; 10; 
         19; 6; 3;];

Costs =[FC;TC];


V_0 = x_0'*Costs;    %  给定的最优解向量 对应的原成本矩阵的花费值


v0 = find(vi == 1);   % 要分别找到vi为 0 或 1 的下标

v1 = find(vi == 0);



model.modelname = 'LB_Inv_UFL';
model.modelsense = 'min';


fi;
% 需要给出 vi 为1的下标  以及 为0 的下标



rik;
 


[h l] = max(reshape(TC,m,n), [], 2);   % 给出每一行的最大值向量m以及下标向量l。  注意1 为每列， 2 为每行。
