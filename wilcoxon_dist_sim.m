function [dist, w] = wilcoxon_dist_sim(N, pH0, ptie, peq0, strateq0)

pfunc = @(Neq0, T) wilcoxon_dist_groups(pH0, Neq0, T, strateq0);

dist = wilcoxon_iterate_groups(pfunc, N, ptie, peq0);
w = (0:0.5:(N*(N+1)/2));  % range of possible values
