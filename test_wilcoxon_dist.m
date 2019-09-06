%% 

N = 10;
pH0 = 0.5;
pH1 = 0.9;
ptie = 0.5;
peq0 = 0.5;
alpha = 0.05;

ptievec = 0:0.1:1.0;
peqvec = 0:0.1:1.0;

%% 

[distbase, wbase] = wilcoxon_dist_simple(N, pH0);
[distwilcox, wwilcox] = wilcoxon_dist_sim(N, pH0, ptie, peq0, 'Wilcoxon');
[distpratt, wpratt] = wilcoxon_dist_sim(N, pH0, ptie, peq0, 'Pratt');
[distmaras, wmaras] = wilcoxon_dist_sim(N, pH0, ptie, peq0, 'Marascuilo');

figure;
subplot(2,1,1);
plot(wbase, distbase, 'k*-');
hold on;
plot(wwilcox, distwilcox, 'o-');
plot(wpratt, distpratt, 'x--');
plot(wmaras, distmaras, '^--');
hold off;
subplot(2,1,2);
stairs(wbase, cumsum(distbase), 'k');
hold on;
stairs(wwilcox, cumsum(distwilcox), '-');
stairs(wpratt, cumsum(distpratt), '--');
stairs(wmaras, cumsum(distmaras), '--');
hold off;

%% 

%%

% wcrit = zeros(length(ptievec), length(peqvec));
% tdx = 1;
% for ptie = ptievec
%   edx = 1;
%   for peq = peqvec
%       [dist, w] = wilcoxon_dist(N, pquart, ptie, peq);
%       wdx = find(cumsum(dist) <= alpha);
%       wcrit(tdx,edx) = w(wdx(end));
%       
%       edx = edx + 1;
%   end
%   tdx = tdx + 1;
% end
% 
% figure;
% imagesc(peqvec, ptievec, wcrit);

