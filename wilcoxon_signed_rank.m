function [p,W] = wilcoxon_signed_rank(x, mu, H0, dim)

%*****************************************************************************
% Copyright (c) 2018      Fiete Winter                                       *
%                         Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock, Germany  *
%                                                                            *
% This file is part of the supplementary material for Fiete Winter's         *
% scientific work and publications                                           *
%                                                                            *
% You can redistribute the material and/or modify it  under the terms of the *
% GNU  General  Public  License as published by the Free Software Foundation *
% , either version 3 of the License,  or (at your option) any later version. *
%                                                                            *
% This Material is distributed in the hope that it will be useful, but       *
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
% or FITNESS FOR A PARTICULAR PURPOSE.                                       *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy of the GNU General Public License along   *
% with this program. If not, see <http://www.gnu.org/licenses/>.             *
%                                                                            *
% http://github.com/fietew/publications           fiete.winter@uni-rostock.de*
%*****************************************************************************

%% Input Parameters
Ndims = ndims(x);
if dim > Ndims
  Ndims = dim;
end
x = permute(x, [dim:Ndims, 1:dim-1]);  % move dim to first dimension
s = size(x);
Nx = s(1);
x = reshape(x, Nx, []);  % squeeze all other dimensions
Nm = size(x,2);

%% Algorithm
x = x-mu;
xabs = abs(x);

pdf = W_dist(Nx);  % distribution of test statistic
Wmax = Nx*(Nx+1)/2;  % naximum value of test statistic

ranks_as = zeros(Nx,1);
ranks_des = zeros(Nx,1);
p = zeros(1,Nm);
W = zeros(1,Nm);
n = 1:Nx;
for idx=1:Nm
    [~, idx_as] = sort(xabs(:,idx),1,'ascend');
    [~, idx_des] = sort(xabs(:,idx),1,'descend');
    
    % get ranks (for ties the average of spanned rank is assigned) 
    ranks_as(idx_as) = n;
    ranks_des(idx_des) = Nx-n+1;
    ranks = 0.5.*ranks_as + 0.5.*ranks_des;
    
    % number of values being greater than mu (0.5 for values equal to mu)
    I= 0.5.*sign(x(:,idx)) + 0.5;
    
    % compute statistic
    W(idx) = sum(ranks.*I);
    
    switch H0
        case '='
            w = ceil(min(W(idx), Wmax - W(idx)));
            p(idx) = 2.*sum(pdf(1:w+1));
        case '<='
            w = ceil(Wmax - W(idx));
            p(idx) = sum(pdf(1:w+1));
        case '>='
            w = ceil(W(idx));
            p(idx) = sum(pdf(1:w+1));
    end
end

%% Output
p = reshape(p, [1, s(2:end)]);  % undo reshape of x
p = permute(p, [Ndims-dim+2:Ndims, 1:Ndims-dim+1]);  % undo permute of x

W = reshape(W, [1, s(2:end)]);  % undo reshape of x
W = permute(W, [Ndims-dim+2:Ndims, 1:Ndims-dim+1]);  % undo permute of x

end

%% Helper Function
function [dist, w] = W_dist(n)
    w = 0:(n*(n+1)/2);  % range of possible values
    dist = zeros(size(w));
    dist(1) = 1;
    for idx=1:n
       dist = 0.5.*dist + 0.5.*[zeros(1,idx),dist(1:end-idx)];
    end
end