function res = wilcoxon_iterate_groups(pfunc, N, ptie, peq0)

res = 0;
for Neq0=0:N  % iterate over possible number of 0's
  pzerocomb = binopdf(Neq0,N,peq0);  % propability of current number of 0's

  % if Neq0 >= N-1, no additional ties can occur
  if Neq0 >= N-1
    res = res + pzerocomb.*pfunc(Neq0, N-Neq0);
  else
    Nt = N-Neq0-1;  % maximum number of ties for current number of 0's
    Ties = zeros(1,Nt);  % array indicating ties between values    
    jj = 1;  % this loop iterates over all possible combinations of ties
    while true
      % if the index of the current parameter is not at its end
      if Ties(jj) <= 1  
        if (jj == Nt)  % if this is the last entry of the array
          % propability of current combination of Ties
          ptiecomb = prod(ptie.^(Ties).*(1-ptie).^(1-Ties));

          if ptiecomb
            % compute groups sizes
            groups = cumsum([1, ~Ties]);
            T = histcounts(groups, groups(end));

            % combine
            res = res + pzerocomb.*ptiecomb.*pfunc(Neq0, T);
          end
          %
          Ties(jj) = Ties(jj) + 1;  % add 1 to the index of the last parameter
        else
          jj = jj + 1;  % if this is not the last entry just move one step
        end
      else
        if jj ~= 1
          Ties(jj) = 0;  % the current index has reached its end, so reset it
          jj = jj - 1;  % move one step backward
          Ties(jj) = Ties(jj) + 1; % ... and add 1 to that index
        else
          % if this is the first index, one can't move backward
          % this means the first index reach its end (you're finished)
          break
        end
      end
    end  % end while   
  end  % end if
end
