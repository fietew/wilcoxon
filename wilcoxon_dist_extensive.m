function [dist, w] = wilcoxon_dist_with_ties(N, pq, ptie)

w = (0:0.5:(N*(N+1)/2)).';  % range of possible values
dist = zeros(length(w),1);

Vars = zeros(2*N-1,1);  % N for Signs, N-1 for Ties
ranks_as = zeros(N,1);
ranks_des = zeros(N,1);

jj = 1;
while true
  % if the index of the current parameter is not at its end
  if Vars(jj) <= 1
    % if this is the last entry
    if (jj == 2*N-1)  % if
      % Actual Code begin
      Signs = Vars(1:N);
      Ties = Vars(N+1:end);
      
      prop = prod(pq.^(Signs).*(1-pq).^(1-Signs)) .* ...
        prod(ptie.^(Ties).*(1-ptie).^(1-Ties));
      
      vals = cumsum([1; ~Ties]);
      [~, idx_as] = sort(vals,1,'ascend');
      [~, idx_des] = sort(vals,1,'descend');
      
      ranks_as(idx_as) = 1:N;
      ranks_des(idx_des) = N-(1:N)+1;
      ranks = 0.5.*ranks_as + 0.5.*ranks_des;
      
      ranksum = Signs.'*ranks;
      idx = ranksum*2+1;
      dist(idx) = dist(idx) + prop;
      % Actual Code end
      %
      Vars(jj) = Vars(jj) + 1;  % add 1 to the index of the last parameter
    else
      jj = jj + 1;  % if this is not the last parameter just move one step
    end
  else
    if jj ~= 1
      Vars(jj) = 0;  % the current index has reached its end, so reset it
      jj = jj - 1;  % move one step backward
      Vars(jj) = Vars(jj) + 1; % ... and add 1 to that index
    else
      % if this is the first index, one can't move backward
      % this means the first index reach its end (you're finished)
      break
    end
  end
end

end

