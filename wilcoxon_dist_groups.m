function [dist, w] = wilcoxon_dist_groups(Neq0, T, pquart, strateq0)

N = Neq0 + sum(T);
dist = zeros(1,N*(N+1)+1);

ranks = Neq0+cumsum(T)-0.5*(T-1);  % ranks for non-zero groups

% handle 0's
if Neq > 0  
  switch strateq0
    case 'Wilcoxon'
      % after Wilcoxon 1945 "Individual Comparisons by Ranking
      % Methods":
      % discard 0's and adjust remaining ranks by subtracting Neq0
      %
      % adjust the ranks
      ranks = ranks - Neq0;
      % The ranksum resulting from the Neq 0's is R = 0
      dist(1) = 1;
    case 'Pratt'
      % after Pratt 1959 "Remarks on Zeros and Ties in the Wilcoxon
      % Signed Rank Procedures":
      % discard 0's and without adjusting the remaining ranks.
      
      % The ranksum resulting from the Neq 0's is R = 0
      dist(1) = 1;
    case 'Marascuilo'
      % after Marascuilo 1977 "Nonparametric and Distribution-free
      % Methods for the Social Sciences":
      % shared rank of the 0's is accounted by one half to the ranksum
      %
      % The first Neq0 values share the rank (1 + 2 + ... + Neq0)/ Neq0
      % = (Neq0+1)/2. Their ranksum is R = Neq*(Neq0+1)/4
      Req0 = Neq0*(Neq0+1)/4;
      dist(2*Req0+1) = 1;
  end
end

% handle remaining ranks for non-zero values
for rdx=1:length(T)
  shift = 2*ranks(rdx);
  for tdx=1:T(rdx)
    dist = pquart*dist + (1-pquart).*[zeros(1,shift),dist(1:end-shift)];
  end
end
