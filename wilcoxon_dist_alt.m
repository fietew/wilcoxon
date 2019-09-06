function [dist, w] = wilcoxon_dist_alt(N, pquart, ptie, peq0, strateq0)

w = (0:0.5:(N*(N+1)/2));  % range of possible values
dist = zeros(size(w));

Ties = zeros(1,N-1);
jj = 1;
while true
  % if the index of the current parameter is not at its end
  if Ties(jj) <= 1
    % if this is the last entry
    if (jj == N-1)  % if
      %%%%%% Actual Code begin
      
      % propability of current combination of Ties
      prop = prod(ptie.^(Ties).*(1-ptie).^(1-Ties));
      
      % compute groups sizes
      groups = cumsum([1, ~Ties]);
      T = histcounts(groups, groups(end));
      
      % Distribution of ranksum R, if there are no 0's
      nonzerodist = wilcoxon_dist_groups(0, T, pquart, strateq0);
      
      % Distribution of ranksum R, if there are 0's
      zerodist = wilcoxon_dist_groups(T(1), T(2:end), pquart, strateq0);
      
      % combine the two distributions
      dist = dist + prop.*((1-peq0).*nonzerodist + peq0.*zerodist);
      %%%%% Actual Code end
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
