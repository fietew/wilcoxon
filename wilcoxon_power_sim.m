function power = wilcoxon_power_sim(N, pH0, pH1, alpha, ptie, peq0, ...
  strateq0, consider_ties_zeros)

pfunc = @(Neq0, T) wilcoxon_power_groups(pH0, pH1, alpha, Neq0, T, ...
  strateq0, consider_ties_zeros);

power = wilcoxon_iterate_groups(pfunc, N, ptie, peq0);
