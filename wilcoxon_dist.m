function [dist, w] = wilcoxon_dist(N, pquart, ptie, peq0, strateq0)

w = (0:0.5:(N*(N+1)/2));  % range of possible values
dist = zeros(size(w));
nonzerodist = zeros(size(w));
zerodist = zeros(size(w));

ranks_as = zeros(1,N);
ranks_des = zeros(1,N);

sdx = find(strcmp(strateq0, {'Wilcoxon', 'Pratt', 'Marascuilo'}));

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
      
      % compute ranks for given ties
      vals = cumsum([1, ~Ties]);
      [~, idx_as] = sort(vals,2,'ascend');
      [~, idx_des] = sort(vals,2,'descend');
      
      ranks_as(idx_as) = 1:N;
      ranks_des(idx_des) = N-(1:N)+1;
      ranks = 0.5.*ranks_as + 0.5.*ranks_des;
      
      %% Distribution of ranksum R, if there are no 0's
      nonzerodist(:) = 0;
      nonzerodist(1) = 1;
      for iidx=2*ranks
        nonzerodist = pquart*nonzerodist ...
          + (1-pquart)*[zeros(1,iidx),nonzerodist(1:end-iidx)];
      end
      
      %% Distribution of ranksum R, if there are 0's
      
      % number of 0's for a given Tie structure
      eqcand = find(~Ties);
      if isempty(eqcand)
        Neq0 = N;
      else
        Neq0 = eqcand(1);
      end
      
      % handle 0's
      zerodist(:) = 0;
      switch sdx
        case 1
          % after Wilcoxon 1945 "Individual Comparisons by Ranking 
          % Methods":
          % discard 0's and adjust remaining ranks by subtracting Neq0
          %
          % adjust the ranks
          ranks = ranks - Neq0;
          % The ranksum resulting from the Neq 0's is R = 0
          zerodist(1) = 1;
        case 2
          % after Pratt 1959 "Remarks on Zeros and Ties in the Wilcoxon
          % Signed Rank Procedures":
          % discard 0's and without adjusting the remaining ranks.
          
          % The ranksum resulting from the Neq 0's is R = 0
          zerodist(1) = 1;
        case 3
          % after Marascuilo 1977 "Nonparametric and Distribution-free
          % Methods for the Social Sciences":
          % shared rank of the 0's is accounted by one half to the ranksum
          %
          % The first Neq0 values share the rank (1 + 2 + ... + Neq0)/ Neq0
          % = (Neq0+1)/2. Their ranksum is R = Neq*(Neq0+1)/4
          Req0 = Neq0*(Neq0+1)/4;
          zerodist(2*Req0+1) = 1;
      end

      % handle remaining ranks for non-zero values
      for iidx=2*ranks(Neq0+1:end)
        zerodist = pquart*zerodist ...
          + (1-pquart)*[zeros(1,iidx),zerodist(1:end-iidx)];
      end
      
      %%
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
