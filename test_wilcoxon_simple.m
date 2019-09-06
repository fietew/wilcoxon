
pquart = 0.5;
alphavec = [0.05,0.025,0.01, 0.005];
Nvec = [6:20, 25:5:50];

wcrit = zeros(length(alphavec), length(Nvec));
adx = 1;
for alpha = alphavec
  ndx = 1;
  for N=Nvec
    [dist, w] = wilcoxon_dist_simple(N, pquart);

    wdx = find(cumsum(dist) <= alpha);
    if isempty(wdx)
        wcrit(adx,ndx) = NaN;
    else
        wcrit(adx,ndx) = w(wdx(end));
    end

    plot(w, cumsum(dist), 'o-');

    hold on;
    ndx = ndx + 1;
  end
  adx = adx + 1;
end

figure;
plot(Nvec, wcrit, 'o--');