function parents = paretotour3(gp, nParents)
% parents: indices (1..popSize) selected for mating (nParents long)
% gp: GPTIPS structure; expects gp.fitness.multiobj to be set (popSize x 3)

popSize = gp.runcontrol.pop_size;
if nargin < 2, nParents = popSize; end

% Prepare matrix of objectives: smaller is better. Order: [fitness, complexity, constraint]
% fitness is RMSE (lower better), complexity currently smaller better, constraint smaller better.
if ~isfield(gp.fitness, 'multiobj')
    error('paretotour3: gp.fitness.multiobj not found. Ensure evalfitness stores multiobj.');
end
P = gp.fitness.multiobj;  % popSize x 3

% Precompute nondominated fronts and ranks
[fronts, rank] = nondominated_sort(P);

parents = zeros(nParents,1);
if isfield(gp.selection,'tournament_size'), tsize = gp.selection.tournament_size; else tsize = 2; end

for s = 1:nParents
    % pick tournament competitors at random
    competitors = randi(popSize, tsize, 1);
    % evaluate which competitor wins
    best = competitors(1);
    for k = 2:tsize
        cand = competitors(k);
        % compare ranks
        if rank(cand) < rank(best)
            best = cand;
        elseif rank(cand) == rank(best)
            % tie: use crowding distance within that front
            fidx = rank(cand);
            front_members = fronts{fidx};
            % get local objective rows and compute cd
            localP = P(front_members,:);
            cd = crowding_distance(localP);
            % map cand and best to local indices
            loc_map = zeros(max(front_members),1); loc_map(front_members) = 1:numel(front_members);
            local_cand_idx = loc_map(cand);
            local_best_idx = loc_map(best);
            % If either cd is Inf (edge), prefer Inf
            cd_c = cd(local_cand_idx);
            cd_b = cd(local_best_idx);
            if cd_c > cd_b
                best = cand;
            end
        end
    end
    parents(s) = best;
end
end
