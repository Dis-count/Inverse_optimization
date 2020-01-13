function [opt,cost_2] = remove0(fi,trans,ind_open)

% This function give best remove [solution,cost] of a given cost (fi,trans)
  [m,n] = size(trans);

  cost_2 = 1e6;
  a = zeros(m,1);
  b = zeros(m*n,1);

  for i = ind_open'  % for (m-p) local circulation

    ini_open = ind_open;

    ini_open(ini_open==i)=[];  % Remove i opening facility

    tran = trans(ini_open,:);  % Extract the selected rows

    [min_tran,tran_index] = min(tran,[],1);  % 1 for column

    remove_cost  = sum(fi(ini_open))+ sum(min_tran);

    if remove_cost < cost_2

      cost_2 = remove_cost;

      ind = ini_open;

      tranindex = tran_index;

    end

  end

    trans_index = [(ind(tranindex)-1)*n]+[1:n]';  % obtain the index of transportation
    a(ind)=1;

    b(trans_index)=1;

    opt = [a;b];  % the suboptimal solution

end
