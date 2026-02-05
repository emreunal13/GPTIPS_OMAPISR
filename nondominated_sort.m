function [fronts, rank] = nondominated_sort(P)
% NONDOMINATED_SORT Simple nondominated sorting (returns cell array fronts and rank vector)
% P is (N x M) where *smaller* values are better.
N = size(P,1);
domCount = zeros(N,1);     % number of solutions that dominate i
S = cell(N,1);             % S{i} stores solutions dominated by i
fronts = {};
rank = inf(N,1);

for i = 1:N
    for j = 1:N
        if i == j, continue; end
        % does i dominate j? i <= j in all objectives AND < in at least one
        pi = P(i,:); pj = P(j,:);
        le = all(pi <= pj);
        lt = any(pi < pj);
        if le && lt
            S{i}(end+1) = j; %#ok<AGROW>
        elseif all(pj <= pi) && any(pj < pi)
            domCount(i) = domCount(i) + 1;
        end
    end
    if domCount(i) == 0
        rank(i) = 1;
    end
end

front = find(rank == 1);
fronts{1} = front;
currFront = 1;

while ~isempty(fronts{currFront})
    Q = []; % next front
    for p = fronts{currFront}
        dominatedList = S{p};
        for q = dominatedList
            domCount(q) = domCount(q) - 1;
            if domCount(q) == 0
                rank(q) = currFront + 1;
                Q(end+1) = q; %#ok<AGROW>
            end
        end
    end
    currFront = currFront + 1;
    fronts{currFront} = unique(Q);
end

% remove final empty cell if present
if isempty(fronts{end})
    fronts(end) = [];
end

end

