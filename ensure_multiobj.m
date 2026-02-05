function gp = ensure_multiobj(gp)
popSize = gp.runcontrol.pop_size;

% pull basics
f = gp.fitness.values(:);
c = gp.fitness.complexity(:);

% constraint score (defaults to 0.5 if missing)
cs = 0.5*ones(popSize,1);
if isfield(gp.fitness,'constraint_values') && ~isempty(gp.fitness.constraint_values)
    n = min(popSize, numel(gp.fitness.constraint_values));
    for i = 1:n
        rv = gp.fitness.constraint_values{i};
        if isstruct(rv)
            if isfield(rv,'Cscore') && isscalar(rv.Cscore) && isfinite(rv.Cscore)
                cs(i) = rv.Cscore;
            else
                acc = [];
                for nm = ["C1","C2","C3","C4"]
                    if isfield(rv, nm) && isscalar(rv.(nm)) && isfinite(rv.(nm))
                        acc(end+1) = rv.(nm); %#ok<AGROW>
                    end
                end
                if ~isempty(acc), cs(i) = mean(acc,'omitnan'); end
            end
        elseif isnumeric(rv) && isscalar(rv) && isfinite(rv)
            cs(i) = rv;
        end
    end
end

% sanitize
f(~isfinite(f))  = 1e6;
c(~isfinite(c))  = 1e6;
cs(~isfinite(cs))= 1;

gp.fitness.multiobj = [f, c, cs];
end

