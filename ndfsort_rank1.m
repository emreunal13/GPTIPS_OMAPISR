function xrank = ndfsort_rank1(x)
%NDFSORT_RANK1 Rank-1 non-dominated solutions (any number of objectives).
%
%   XRANK = NDFSORT_RANK1(X) takes an (N x M) objective array X where
%   *smaller values are better* in every column, performs non-dominated
%   sorting, and returns a column vector XRANK of length N with
%
%       XRANK(i) = 1  if solution i is on the (global) Pareto front
%       XRANK(i) = 0  otherwise.
%
%   Historically this was restricted to 2 objectives; now it delegates
%   to the generic N-objective NONDOMINATED_SORT function.
%
%   See also NONDOMINATED_SORT.

    N = size(x,1);
    xrank = zeros(N,1);

    if N == 0
        return;
    end

    [fronts, ~] = nondominated_sort(x);  % expects "lower is better"

    if ~isempty(fronts) && ~isempty(fronts{1})
        xrank(fronts{1}) = 1;
    end
end
