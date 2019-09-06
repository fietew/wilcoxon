function res = wilcoxon_iterate_groups(pfunc, N, ptie, peq0)

res = 0;
Ties = zeros(1,N-1);
jj = 1;
while true
  % if the index of the current parameter is not at its end
  if Ties(jj) <= 1
    % if this is the last entry
    if (jj == N-1)  % if
      % propability of current combination of Ties
      prop = prod(ptie.^(Ties).*(1-ptie).^(1-Ties));
      
      if prop
        % compute groups sizes
        groups = cumsum([1, ~Ties]);
        T = histcounts(groups, groups(end));
        
        % combine
        res = res + prop* ...
          ((1-peq0)*pfunc(0, T) + peq0*pfunc(T(1), T(2:end)));
      end
      %
      Ties(jj) = Ties(jj) + 1;  % add 1 to the index of the last parameter
    else
      jj = jj + 1;  % if this is not the last parameter just move one step
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
end
