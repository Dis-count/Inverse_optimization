% This function is used to get excel to show the specific data result.
function Ran_UFL(m,n,k)   
% (m,n) 问题规模  k 为调用函数次数 默认 k=50;

theresult=zeros(k,2);

for j = 2:10 

    for i=1:k
    
        [a,b] =Random_UFL(m,n,j);  % 重复调用
        theresult(i,:) = [a,b];
    end


filename = ['F:\Program Files\Matlab files\',num2str(m),'by',num2str(n),'-',num2str(j),'.xlsx'];

xlswrite(filename,theresult);

end
