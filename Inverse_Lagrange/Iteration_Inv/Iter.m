function A0 = Iter(iter_max)

% This function is used to iterate to obtain the value of the corresponding
A = [6 4;-1 1];
c = [5;4];
b = [24;1];
x0 = [2;2.5];

y = [2/3;4];

i = 1;

while (abs(y(1)-5/6)>1e-5)

    xy = Master(A,b,c,y,x0);   % 更新 y

    At = reshape(A',4,1);

    l_A = length(At);

    A0 = xy(1:l_A) - xy(l_A+1:2*l_A) + At;

    A1 = reshape(A0,2,2);

    A0 = A1';

    y0 = Sub_Max(A0,b,c);  % 更新 A0

    y  = (c' * x0)/(b'*y0)*y0;

    i = i + 1;

    if i > iter_max

      break

    end

end

end
