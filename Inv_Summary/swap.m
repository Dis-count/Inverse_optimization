function [opt,cost_1] = swap(fi,trans,ind_open,ind_unopen)
% This function give best swap [solution,cost] of a given cost (fi,trans)

  [m,n] = size(trans);

  cost_1 = 1e6;
  a = zeros(m,1);
  b = zeros(m*n,1);

  for i = ind_unopen'  % for p local circulation

    for j = ind_open'  % for (m-p) local circulation

      ini_open = ind_open;

      ini_open = [ini_open;i];  % Add i opening facility

      ini_open(ini_open==j)=[];   % Remove j opening facility

      ini_open = sort(ini_open);

      tran = trans(ini_open,:);  % Extract the selected rows

      [min_tran,tran_index] = min(tran,[],1);  % 1 for column

      swap_cost  = sum(fi(ini_open))+ sum(min_tran);

      if swap_cost < cost_1

        cost_1 = swap_cost;

        ind = ini_open;  % The optimal opening facility

        tranindex = tran_index;

      end
    end
  end

    trans_index = [(ind(tranindex)-1)*n]+[1:n]'; % Attention the order
     % obtain the index of transportation
    a(ind)=1;

    b(trans_index)=1;

    opt = [a;b];  % the suboptimal solution

end
