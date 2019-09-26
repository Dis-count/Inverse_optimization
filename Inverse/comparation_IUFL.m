m = 6;

n = 6;
k = 10;

res = zeros(1,k);

for j = 1:k

    vi = zeros(m,1);

    cee = [];

    for i = 1: n
        ce = [ones(1,1),zeros(1,m-1)];

        ce = ce(randperm(m));

        cee = [cee;ce];

    end

    c = cee';

    for i = 1: m

      a = logical(sum(c(i,:)));

      vi(i) = a + 0;

    end

    uik = reshape(cee, m*n, 1); % 直接 按列 reshape ��?以不用再转置

    v_UFL = round((-0.5+rand(1,1))*10); % (-1,1)取�?? ��?优�?�的漂移��?

    FC = round(10 + rand(m,1)*50);

    TC = round(5 + rand(m*n,1)*20);

    [opt1,opt2] = UFL(FC, TC);

    LB_res = LB_IUFL(v_UFL + opt1, vi, uik, FC, TC);

    L_res = L_IUFL(v_UFL + opt1, vi, uik, FC, TC);

    % Inv_res = Inv_UFL3(v_UFL+opt1, vi, uik, FC, TC);

    res(j) = (- L_res + LB_res)/ L_res;

end

res
res1 = mean(abs(res))
