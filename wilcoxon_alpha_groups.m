function alpha_actual = wilcoxon_alpha_groups(pH0, alpha, Neq0, T, strateq0, ...
  consider_ties_zeros)

N = Neq0 + sum(T);

% actual distribution (H0 true) of ranksum, if 0's and ties are present
dist_actual = wilcoxon_dist_groups(pH0, Neq0, T, strateq0);

% distribution (H0 true) of ranksum assumed during the test
if consider_ties_zeros
  dist_assumed = dist_actual;
else
  dist_assumed = wilcoxon_dist_groups(pH0, 0, ones(1,N), strateq0);
end

% index of critical value wcrit for rank sum (dist_actual(w <= wcrit) = alpha)
wdx = find(cumsum(dist_assumed) <= alpha, 1, 'last');
if isempty(wdx)
  alpha_actual = 0;
else
  P_actual = cumsum(dist_actual);
  alpha_actual = P_actual(wdx);
end
