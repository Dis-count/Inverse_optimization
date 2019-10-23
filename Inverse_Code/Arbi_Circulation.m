% varRange中数组取值对顺序敏感。
% 循环变量遍历集若是整体(而不是标量数字，或者单个字符等）构成，那么需要用花括号括起来。
varRange = {1:10, {[3 2], -6}, 'abc',{'China','UK'}}; % 循环变量集合，注意加花括号和不加花括号最后结果的区别
nloop = length(varRange);   % for循环层数，亦即循环变量个数
var = cell(nloop,1);    % 初始化循环变量
maxIter = cellfun(@length,varRange(:));
subs = arrayfun(@(x)1:x,maxIter,'un',0);
iter = cell(nloop,1);
[iter{end:-1:1}] = ndgrid(subs{end:-1:1});
temp = cellfun(@(x)x(:),iter,'un',0);
iter = num2cell(cat(2,temp{:}));    % 循环变量索引
for k = 1:size(iter,1)  % 一次循环替代嵌套
    var = cellfun(@(x,y)x(y),varRange,iter(k,:),'un',0);
    % do something
    dispvar(k,nloop,var)
end



% 上面的dispvar是一个用于提示循环变量取值的函数

function dispvar(k,nloop,var)
% k - 当前循环计数
% nloop - 嵌套循环总层数
% var - 当前循环变量取值组合
% 未做详细的参数检查，在非常规输入情形下可能会有bug

fprintf(['循环变量遍历组合%' num2str(fix(log10(k))+1) 'd:??('],k)
for s = 1:nloop
    if ~iscell(var{s})
        if isnumeric(var{s})
            fprintf('%g??',var{s})
        elseif ischar(var{s})
            fprintf('%c??',var{s})
        end
    else
        if isnumeric(var{s}{:}) && ~isscalar(var{s}{:})
            fprintf('[%0.f-by-%0.f %s] ',size(var{s}{:}),class(var{s}{:}));
        elseif isnumeric(var{s}{:}) && isscalar(var{s}{:})
            fprintf('%g ',var{s}{:});
        elseif ischar(var{s}{:})
            fprintf('%s??',var{s}{:})
        end
    end
end
fprintf(')\n')