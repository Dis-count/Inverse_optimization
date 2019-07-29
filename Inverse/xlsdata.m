count = 1;   % This function is used to integrate all the data from excels.

result = zeros(63,3);

for i = 10:10:70
    
    for j =  2:10

    filename = ['F:\Program Files\Matlab files\',num2str(i),'by',num2str(2*i),'-',num2str(j),'.xlsx'];
    
    sheet = 1;
    xlRange = 'A51:C51';

    result(count,:) = xlsread(filename,sheet,xlRange);
    count = count+1;
    end
    
end

filename1 = ['F:\Program Files\Matlab files\theresult.xlsx'];


xlswrite(filename1,result');