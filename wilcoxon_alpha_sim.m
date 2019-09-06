function alpha_actual = wilcoxon_alpha_sim(N, pH0, alpha, ptie, peq0, ...
  strateq0, consider_ties_zeros)

pfunc = @(Neq0, T) wilcoxon_alpha_groups(pH0, alpha, Neq0, T, ...
  strateq0, consider_ties_zeros);

alpha_actual = wilcoxon_iterate_groups(pfunc, N, ptie, peq0);
