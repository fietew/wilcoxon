%%

Nvec = 5:1:12;
pH0 = 0.5;
pH1 = 0.95;
ptie = 0.1;
peq0 = 0.0;
alpha = 0.025;

ptievec = 0:0.1:1.0;
peqvec = 0:0.1:1.0;

%%

figure;
legendstring = {};
for strateq = {'Wilcoxon', 'Pratt', 'Marascuilo'}
  for consider_ties_zeros = 0:1
    power = zeros(size(Nvec));
    alpha_actual = zeros(size(Nvec));
    ndx = 1;
    for N=Nvec
      power(ndx) = wilcoxon_power_sim(N, pH0, pH1, alpha, ptie, peq0, strateq{1}, consider_ties_zeros);
      alpha_actual(ndx) = wilcoxon_alpha_sim(N, pH0, alpha, ptie, peq0, strateq{1}, consider_ties_zeros);
      
      ndx = ndx + 1;
    end
    subplot(3,1,1);
    hold on;
    plot(Nvec, power, 'o-');
    
    subplot(3,1,2);
    hold on;
    plot(Nvec, alpha_actual, 'o-');
    
    subplot(3,1,3);
    hold on;
    plot(alpha_actual, power, 'o-');
    
    legendstring = [legendstring, ...
      sprintf('%s consider:%d', strateq{1}, consider_ties_zeros)];
  end
end
hold off;

subplot(3,1,3);
legend(legendstring, 'Location', 'southeast');
xlim([0,1]);
ylim([0,1]);

