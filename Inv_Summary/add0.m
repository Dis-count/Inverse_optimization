function [opt,cost_0]= add0(fi,trans,ind_open,ind_unopen)
% This function give best add [solution,cost] of a given cost (fi,trans)

  [m,n] = size(trans);

  cost_0 = 1e6;
  a = zeros(m,1);
  b = zeros(m*n,1);

  for j = ind_unopen'  % for p local circulation

    ini_open = ind_open;

    ini_open = [ini_open;j];  % Add the new i opening facility

    ini_open = sort(ini_open);

    tran = trans(ini_open,:);  % Extract the selected rows

    [min_tran,tran_index] = min(tran,[],1);  % 1 for column

    add_cost  = sum(fi(ini_open))+ sum(min_tran);

    if add_cost < cost_0

      cost_0 = add_cost;  % Save the optimal cost

      ind = ini_open;

      tranindex = tran_index;  % Save the optimal index

    end

  end

    trans_index = [(ind(tranindex)-1)*n]+[1:n]';  % obtain the index of transportation
    a(ind)=1;

    b(trans_index)=1;

    opt = [a;b];  % the suboptimal solution

end
