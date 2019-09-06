function [dist, w] = wilcoxon_dist_simple(N, pquart)

    w = (0:(N*(N+1)/2));  % range of possible values
    dist = zeros(size(w));
    dist(1) = 1;
    for idx=1:N
       dist = pquart*dist + (1-pquart)*[zeros(1,idx), dist(1:end-idx)];
    end
end

