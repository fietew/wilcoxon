function power = wilcoxon_power_groups(pH0, pH1, alpha, Neq0, T, ...
  strateq0, consider_ties_zeros)

N = Neq0 + sum(T);

% distribution (H0 true) of ranksum assumed during the test to determine 
% the critical value 
if consider_ties_zeros
  distH0 = wilcoxon_dist_groups(pH0, Neq0, T, strateq0);
else
  distH0 = wilcoxon_dist_groups(pH0, 0, ones(1,N), strateq0);
end
% index of critical value wcrit for rank sum (distH0(w <= wcrit) = alpha)
wdx = find(cumsum(distH0) <= alpha, 1, 'last');
if isempty(wdx)
    power = 0;
    return
end
% distribution of ranksum, if H1 is true
distH1 = wilcoxon_dist_groups(pH1, 0, [Neq0,T], strateq0);
PH1 = cumsum(distH1);
power = PH1(wdx);