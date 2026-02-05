function cd = crowding_distance(P)
% CROWDING_DISTANCE standard NSGA-II crowding distance for rows of P.
% P: nPoints x nObjectives, larger cd means less crowded
n = size(P,1);
m = size(P,2);
cd = zeros(n,1);
if n == 0
    return;
elseif n == 1
    cd(1) = Inf;
    return;
end

for j = 1:m
    [sortedVals, idx] = sort(P(:,j),'ascend');
    cd(idx(1)) = Inf;
    cd(idx(end)) = Inf;
    vmax = sortedVals(end);
    vmin = sortedVals(1);
    if vmax == vmin
        continue;
    end
    % normalized difference for interior points
    for k = 2:(n-1)
        cd(idx(k)) = cd(idx(k)) + (sortedVals(k+1) - sortedVals(k-1)) / (vmax - vmin);
    end
end

end


