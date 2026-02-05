function [fitness,gp,theta,ypredtrain,fitnessTest,ypredtest,pvals,...
          r2train,r2test,r2val,geneOutputs,geneOutputsTest,geneOutputsVal] = ...
          regressmulti_fitfun(evalstr_in,gp)
%REGRESSMULTI_FITFUN Physics-aware, full-constant Nelder–Mead optimizer.
%
%   This replaces the default GPTips multigene fitness function.
%
%   - Optimizes ALL constants at once:
%       * bias + gene weights
%       * all numeric literals inside trees (ERCs, exponents, etc.)
%   - Inner loop objective = RMSE + lambda_phys * fast physics penalty
%   - After NM:
%       * stores theta in gp.fitness.returnvalues{ci}
%       * stores full parameter vector in gp.fitness.param_values{ci}
%       * (hook) optional Lamarckian update of gp.pop{ci}
%       * computes robust C1..C4 and saves in gp.fitness.constraint_values{ci}
%
% Signature kept compatible with original for post-run use.

% ---- defaults in case of early exit ----
theta=[]; ypredtrain=[]; fitnessTest=[]; ypredtest=[]; pvals=[];
r2train=[]; r2test=[]; r2val=[];
geneOutputs=[]; geneOutputsTest=[]; geneOutputsVal=[];

ci = gp.state.current_individual;
y  = gp.userdata.ytrain;
numData = gp.userdata.numytrain;

% Safety
if isempty(y) || numData == 0
    fitness = Inf;
    return;
end

% 0) Build parametric templates: x1,x2-forms and training / X(:,1:2) forms
%   evalstr_in: cellstr with 'x1','x2',... from tree2evalstr
evalstr_raw = evalstr_in;

% Encode all numeric literals as "C(k)" and collect their initial values
[geneTemplates_raw, const0] = encode_numeric_constants_all(evalstr_raw);

numGenes = numel(geneTemplates_raw);
nTheta   = numGenes + 1;        % bias + one weight per gene
nConst   = numel(const0);
D        = numData;

% Build training templates: x1->gp.userdata.xtrain(:,1)
patX = 'x(\d+)';
geneTemplatesTrain = cell(size(geneTemplates_raw));
for g = 1:numGenes
    geneTemplatesTrain{g} = regexprep(geneTemplates_raw{g},patX,...
        'gp.userdata.xtrain(:,$1)');
end

% Build generic X-templates: x1->X(:,1) (for Re,e/D probes)
geneTemplatesX = cell(size(geneTemplates_raw));
for g = 1:numGenes
    geneTemplatesX{g} = regexprep(geneTemplates_raw{g},patX,'X(:,$1)');
end

% 1) Initial parameter vector p0  = [theta_0 ; const0]
%    theta_0 via robust ridge regression on current constants.
phys = get_phys_params_safe(gp);

% If we already have a param vector for this individual, try warm start
p0 = [];
if isfield(gp.fitness,'param_values') && ...
   numel(gp.fitness.param_values) >= ci && ...
   ~isempty(gp.fitness.param_values{ci})

    p_prev = gp.fitness.param_values{ci};
    if numel(p_prev) == nTheta + nConst
        p0 = p_prev(:);
    end
end

if isempty(p0)
    % Build gene output matrix G with current constants const0
    C = const0(:).';
    G = ones(D, numGenes+1);
    for g = 1:numGenes
        G(:,g+1) = eval(geneTemplatesTrain{g});  % now without try/catch
    end

    % Clean up non-finite gene outputs
    badCols = any(~isfinite(G), 1);
    if any(badCols)
        if gp.runcontrol.verbose > 1
            warning('Non-finite gene outputs in individual %d. Zeroing %d columns.', ...
                ci, sum(badCols));
        end
        G(:, badCols) = 0;   % those genes effectively get zero weight
    end


    % Robust ridge LS for theta_0
    if isfield(gp.userdata,'bootSample') && gp.userdata.bootSample
        sampleInds = bootsample(G,gp.userdata.bootSampleSize);
        X = G(sampleInds,:);
        yt = y(sampleInds);
    else
        X = G; yt = y;
    end

    XTX = X.' * X;
    Xty = X.' * yt;
    dim = size(XTX,1);
    avg_energy = trace(XTX)/max(dim,1);
    if avg_energy == 0, avg_energy = 1; end
    lambda_scaled = 1e-3 * avg_energy;

    try
        theta_0 = (XTX + lambda_scaled*eye(dim)) \ Xty;
    catch ME
        %warning('regressmulti_fitfun:ridgeCrash', ...
            %'Ridge LS crashed for individual %d: %s. Using pinv(G)*y.', ...
            %ci, ME.message);
        theta_0 = pinv(G) * y;
    end

    if any(~isfinite(theta_0))
        %warning('Non-finite theta_0 for individual %d. Marking fitness = Inf.', ci);
        fitness = Inf;
        gp.fitness.returnvalues{ci} = [];
        gp = ensure_multiobj_field(gp);
        gp.fitness.multiobj(ci,:) = [Inf, gp.fitness.complexity(ci), 1];
        gp.fitness.constraint_values{ci} = ...
            struct('C1',1,'C2',1,'C3',1,'C4',1,'Cscore',1);
        return;
    end

    % Pack initial vector
    theta_0 = theta_0(:);
    if numel(theta_0) ~= nTheta
        % pad or truncate to match number of genes
        theta_0 = [theta_0(:); zeros(nTheta-numel(theta_0),1)];
        theta_0 = theta_0(1:nTheta);
    end
    p0 = [theta_0; const0(:)];
end


% 2) Nelder–Mead: optimize ALL params p = [theta ; C]
options = optimset('Display','off','MaxIter',800,'MaxFunEvals',2000,'TolFun',1e-10);
lambda_phys = 0;

lossFun = @(p) physics_loss_all_param( ...
    p, gp, geneTemplatesTrain, geneTemplatesX, y, lambda_phys, phys);

try
    p_best = fminsearch(lossFun, p0, options);
catch ME
    warning('regressmulti_fitfun:fminsearchCrash', ...
        'fminsearch crashed for individual %d: %s. Using p0 as p_best.', ...
        ci, ME.message);
    p_best = p0;
end

theta = p_best(1:nTheta);
Cbest = p_best(nTheta+1:end);


% 3) Final train predictions & fitness
C = Cbest(:).';
geneOutputs = ones(D, numGenes+1);

% Any error here is a structural bug -> let it throw
for g = 1:numGenes
    geneOutputs(:,g+1) = eval(geneTemplatesTrain{g});  % uses C, xtrain
end

% Clean nasty columns *before* computing ypred
badCols = any(~isfinite(geneOutputs), 1);
if any(badCols)
    %warning('regressmulti_fitfun:nonFiniteGeneOutputsFinal', ...
        %'Non-finite gene outputs after NM for individual %d. Zeroing %d columns.', ...
        %ci, sum(badCols));
    geneOutputs(:, badCols) = 0;
end

ypredtrain = geneOutputs * theta;
err = (y - ypredtrain) ./ y;

if any(~isfinite(ypredtrain)) || any(~isfinite(err))
    %warning('regressmulti_fitfun:nonFiniteTrainPred', ...
        %'Non-finite ypred/err for individual %d. Setting fitness = 1.', ci);
    fitness = 1;
else
    if max(abs(err)) < 1
        %RMSE
        fitness = sqrt(mean(abs(err.^2)));
    else
        %warning('regressmulti_fitfun:largeRelError', ...
            %'Relative error exceeded 1 for individual %d. Setting fitness = 1.', ci);
        fitness = 1;
    end
end


% 4) Physics constraints C1..C4 and aggregate constraint score
eD_list = [1.87569231e-06, 0.00100603621730382, 0.00202429149797571,...
           0.00404203718674212, 0.00821692686935086, 0.0164338537387017,...
           0.0331674958540630];

c1 = 1; c2 = 1; c3 = 1; c4 = 1;
fitness_gate = 0.7;

if isfinite(fitness) && (fitness <= fitness_gate)
    % Single model handle: f_handle(Re,eD) = friction factor
    f_handle = assemble_model_handle_param(geneTemplatesX, p_best);

    % All four constraints use the same interface (handle + phys)
    c1 = compute_C1_penalty_param(f_handle, phys, eD_list);
    c2 = compute_C2_penalty_param(f_handle, phys, eD_list);
    c3 = compute_C3_penalty_param(f_handle, phys);
    c4 = compute_C4_penalty_param(f_handle, phys);
end

% Aggregate with L2-type mean, clamp to [0,1]
p_agg  = 2;
Cvec   = [c1,c2,c3,c4];

c_score = (mean(Cvec.^p_agg,'omitnan')).^(1/p_agg);
constraint_score = min(max(real(c_score),0),1);

% 5) Store results + Lamarckian hooks
ci = gp.state.current_individual;

gp.fitness.returnvalues{ci}      = theta;   % bias + gene weights
gp.fitness.param_values{ci}      = p_best;  % [theta ; C]
gp.fitness.constraint_values{ci} = struct('C1',c1,'C2',c2,'C3',c3,'C4',c4,...
                                          'Cscore',constraint_score);

% Write optimized constants back into the genotype
% const_new (Cbest) must be in the SAME global order as const0
gp.pop{ci} = update_individual_constants_in_pop(gp.pop{ci}, Cbest);

% 6) Post-run stats / plotting
fitnessTest = []; ypredtest = []; r2train = []; r2test = []; r2val = []; pvals = [];
geneOutputsTest = []; geneOutputsVal = [];

if gp.state.run_completed

    r2train = 1 - sum( (y-ypredtrain).^2 ) / sum( (y-mean(y)).^2 );
    plotValidation = 0; plotTest = 0;

    % validation set 
    if isfield(gp.userdata,'xval') && isfield(gp.userdata,'yval') && ...
       ~isempty(gp.userdata.xval) && ~isempty(gp.userdata.yval)

        plotValidation = 1;
        numDataVal = length(gp.userdata.yval);
        geneOutputsVal = ones(numDataVal,numGenes+1);

        % reuse raw templates: x1->gp.userdata.xval(:,1)
        geneTemplatesVal = cell(size(geneTemplates_raw));
        for g = 1:numGenes
            s = regexprep(geneTemplates_raw{g},'x(\d+)','gp.userdata.xval(:,$1)');
            geneTemplatesVal{g} = s;
            C = Cbest(:).';
            geneOutputsVal(:,g+1) = eval(s);
        end

        ypredval = geneOutputsVal*theta;
        fitness_val = sqrt(mean((gp.userdata.yval - ypredval).^2));
        r2val = 1 - sum( (gp.userdata.yval - ypredval).^2 ) / ...
                    sum( (gp.userdata.yval - mean(gp.userdata.yval)).^2 );
    end

    % test set 
    if isfield(gp.userdata,'xtest') && isfield(gp.userdata,'ytest') && ...
       ~isempty(gp.userdata.xtest) && ~isempty(gp.userdata.ytest)

        plotTest = 1;
        numDataTest = length(gp.userdata.ytest);
        geneOutputsTest = ones(numDataTest,numGenes+1);

        geneTemplatesTest = cell(size(geneTemplates_raw));
        for g = 1:numGenes
            s = regexprep(geneTemplates_raw{g},'x(\d+)','gp.userdata.xtest(:,$1)');
            geneTemplatesTest{g} = s;
            C = Cbest(:).';
            geneOutputsTest(:,g+1) = eval(s);
        end

        ypredtest = geneOutputsTest * theta;
        fitnessTest = sqrt(mean((gp.userdata.ytest - ypredtest).^2));
        r2test = 1 - sum( (gp.userdata.ytest - ypredtest).^2 ) / ...
                      sum( (gp.userdata.ytest - mean(gp.userdata.ytest)).^2 );
    end

    % stats toolbox p-values for genes 
    if gp.userdata.stats && gp.info.toolbox.stats
        wstate = warning; warning off;
        stats = regstats(y,geneOutputs(:,2:end));
        warning(wstate);
        pvals = stats.tstat.pval;
    end

    % plots
    if gp.userdata.showgraphs
        plot_regressmulti_graphs(gp, y, ypredtrain, fitness, r2train, ...
                                 plotTest, ypredtest, fitnessTest, r2test, ...
                                 plotValidation, ypredval, fitness_val, r2val, ...
                                 numGenes, geneOutputs, pvals);
    end
end

end % main function


%% Utility: ensure multiobj field 
function gp = ensure_multiobj_field(gp)
    popSize = gp.runcontrol.pop_size;
    if ~isfield(gp.fitness,'multiobj') || isempty(gp.fitness.multiobj) ...
            || size(gp.fitness.multiobj,1) ~= popSize
        gp.fitness.multiobj = nan(popSize,3);
    end
end

%% Encode numeric constants as C(k) 
function [templates, const_vals] = encode_numeric_constants_all(evalstr_raw)
% evalstr_raw: cell of strings using x1,x2,...
% Output:
%   templates : same, but all numeric literals -> C(k)
%   const_vals: vector of those numeric literals in order of appearance

    templates   = evalstr_raw;
    const_vals  = [];
    global_idx  = 0;

    % Regex for standalone numeric literal
    numPattern = '(?<![A-Za-z0-9_\.])([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';

    for g = 1:numel(evalstr_raw)
        s = evalstr_raw{g};
        out = '';
        lastIdx = 1;
        [startIdx, endIdx, ~, matches] = regexp(s, numPattern, 'start','end','tokenExtents','tokens');

        if isempty(startIdx)
            templates{g} = s;
            continue;
        end

        for m = 1:numel(startIdx)
            % literal numeric text
            lit = matches{m}{1};
            global_idx = global_idx + 1;
            const_vals(global_idx,1) = str2double(lit);

            out = [out, s(lastIdx:startIdx(m)-1), 'C(', num2str(global_idx), ')'];
            lastIdx = endIdx(m)+1;
        end
        out = [out, s(lastIdx:end)]; %#ok<AGROW>
        templates{g} = out;
    end
end

%% Physics params from gp 
function phys = get_phys_params_safe(gp)
    phys.D      = 0.05;
    phys.L      = 1.0;
    phys.mu_ref = 1e-3;
    phys.rho    = 1000;

    if isfield(gp,'userdata') && isstruct(gp.userdata)
        if isfield(gp.userdata,'D')      && ~isempty(gp.userdata.D),      phys.D      = gp.userdata.D; end
        if isfield(gp.userdata,'L')      && ~isempty(gp.userdata.L),      phys.L      = gp.userdata.L; end
        if isfield(gp.userdata,'mu_ref') && ~isempty(gp.userdata.mu_ref), phys.mu_ref = gp.userdata.mu_ref; end
        if isfield(gp.userdata,'rho')    && ~isempty(gp.userdata.rho),    phys.rho    = gp.userdata.rho; end
    end
    if isfield(gp,'user') && isstruct(gp.user)
        if isfield(gp.user,'D')      && ~isempty(gp.user.D),      phys.D      = gp.user.D; end
        if isfield(gp.user,'L')      && ~isempty(gp.user.L),      phys.L      = gp.user.L; end
        if isfield(gp.user,'mu_ref') && ~isempty(gp.user.mu_ref), phys.mu_ref = gp.user.mu_ref; end
        if isfield(gp.user,'rho')    && ~isempty(gp.user.rho),    phys.rho    = gp.user.rho; end
    end
end

function cost = physics_loss_all_param(p, gp, geneTemplatesTrain, geneTemplatesX, y, lambda, phys)
    % Inner objective for NM: RMSE + lambda * fast-physics penalty.
    % Now has access to *gp* so geneTemplatesTrain that reference
    % gp.userdata.xtrain actually work.

    numGenes = numel(geneTemplatesTrain);
    nTheta   = numGenes + 1;

    theta = p(1:nTheta);
    C     = p(nTheta+1:end);
    C     = C(:).';  % row for eval(C(k))

    % ---- build gene output matrix on training data ----
    D = gp.userdata.numytrain;
    G = ones(D, numGenes+1);

    for g = 1:numGenes
        try
            G(:,g+1) = eval(geneTemplatesTrain{g});
        catch ME
            error('physics_loss_all_param: invalid gene template (ind %d, gene %d): "%s"\n%s', ...
                gp.state.current_individual, g, geneTemplatesTrain{g}, ME.message);
        end
    end

    if any(~isfinite(G(:)))
        cost = 1e6;   
        return;
    end

    ypred = G * theta;
    if any(~isfinite(ypred))
        cost = 1e6;
        return;
    end

    % if we want fully real ypred:
    % ypred = real(ypred);

    err  = (y - ypred) ./ y;
    rmse = sqrt(mean(abs(err).^2));   % always real 

    % ---- fast physics penalty ----
    %pen_fast = compute_fast_physics_penalty_param(geneTemplatesX, p, phys);
    pen_fast = 0;

    cost = rmse + lambda * pen_fast;
end

%% Generic model eval on X = [Re,eD]
function f = eval_model_param(geneTemplatesX, theta, C, X)
    numGenes = numel(geneTemplatesX);
    C = C(:).';
    G = ones(size(X,1), numGenes+1);
    for g = 1:numGenes
        
        G(:,g+1) = eval(geneTemplatesX{g});  % uses X, C
    end
    f = G * theta;
end


function pen = compute_C1_penalty_param(f_handle, phys, eD_list)
% C1 in [0,1] (0 = best).
%   For each e/D:
%     - Zone A: 1e4 <= Re <= 1e6  →  1 <= n <= 2.4
%     - Zone B: 1e6 <= Re <= 1e9  →  n ≈ 2
%
%   n = d log(ΔP) / d log(V) with ΔP = f * ρ V^2 L / (2D)

    rho    = phys.rho;
    mu_ref = phys.mu_ref;
    D      = phys.D;
    L      = phys.L;

    % Re grid 
    Re_grid = logspace(4, 9, 160).';     % 1e4 .. 1e9
    V_grid  = (Re_grid .* mu_ref) ./ (rho * D);

    ne = numel(eD_list);
    per_eD = nan(ne,1);

    movmean3 = @(v) conv(v,[1 1 1]/3,'same');

    % zone definitions
    ReA_min = 1e4; ReA_max = 1e6;
    ReB_min = 1e6; ReB_max = 1e9;

    nA_min  = 1.0;
    nA_max  = 2.4;
    nB_ref  = 2.0;
    tolB    = 0.001;   % how close to 2 we want n in zone B

    for ii = 1:ne
        eD = eD_list(ii);

        % f(Re,e/D)
        f = f_handle(Re_grid, eD*ones(size(Re_grid)));

        ok = isfinite(f) & isreal(f) & (f > 0);
        if nnz(ok) < 30
            per_eD(ii) = 1;   % not enough valid points
            continue;
        end

        Re_ok = Re_grid(ok);
        V_ok  = V_grid(ok);
        f_ok  = f(ok);

        % deltaP vs V
        DP_ok = f_ok .* (rho .* (V_ok.^2) .* L) ./ (2 .* D);
        ok2 = isfinite(DP_ok) & isreal(DP_ok) & (DP_ok > 0);
        if nnz(ok2) < 30
            per_eD(ii) = 1;
            continue;
        end

        Re_ok = Re_ok(ok2);
        V_ok  = V_ok(ok2);
        DP_ok = DP_ok(ok2);

        % local n(V)
        if numel(DP_ok) >= 5
            DP_sm = movmean3(DP_ok);
            V_sm  = movmean3(V_ok);
        else
            DP_sm = DP_ok;
            V_sm  = V_ok;
        end

        lnDP = log(DP_sm);
        lnV  = log(V_sm);
        nloc = gradient(lnDP, lnV);   % n(V)

        % Zone A: 1 <= n <= 2.4  (mixed / transitional) 
        maskA = (Re_ok >= ReA_min) & (Re_ok <= ReA_max);
        if nnz(maskA) < 6
            penA = 1;   % if we can't assess, treat as bad
        else
            nA = nloc(maskA);

            dev_low  = max(0, nA_min - nA);   % violations below 1
            dev_high = max(0, nA    - nA_max);% violations above 2.4

            % normalise by band width (2.4 - 1 = 1.4) and clip
            devA = mean(dev_low + dev_high, 'omitnan') / (nA_max - nA_min + eps);
            penA = min(max(real(devA),0),1);
        end

        % Zone B: n ≈ 2 (fully rough tail) 
        maskB = (Re_ok >= ReB_min) & (Re_ok <= ReB_max);
        if nnz(maskB) < 6
            penB = 1;
        else
            nB = nloc(maskB);
            devB_raw = mean( max(0, abs(nB - nB_ref) - tolB), 'omitnan' );
            % scale so ~0.5 away -> near penalty 1
            devB = devB_raw / (0.5 + eps);
            penB = min(max(real(devB),0),1);
        end

        % combine zones for this e/D
        wA = 0.3; wB = 0.7;
        per_eD(ii) = min(1, wA*penA + wB*penB);
    end

    % overall C1 = worst e/D
    C1_core = max(per_eD, [], 'omitnan');
    if ~isfinite(C1_core), C1_core = 1; end

    pen = min(max(real(C1_core),0),1);
end


function pen = compute_C2_penalty_param(f_handle, phys, eD_list)
% C2 in [0,1] (0 = best).
% deltaP(e) should be weakly non-decreasing at fixed Re, with small dips allowed.

rho    = phys.rho;
mu_ref = phys.mu_ref;
D      = phys.D;
L      = phys.L;

Re_probe = logspace(log10(3e3), log10(1e8), 12);
egrid    = logspace(log10(min(eD_list)*D*0.5), ...
                    log10(max(eD_list)*D*1.5), 360).';

smooth_w = 5;
s_allow  = 0.05;   % allow s(e) >= -0.05 roughly

viol_all = nan(numel(Re_probe),1);

for k = 1:numel(Re_probe)
    Re0 = Re_probe(k);

    % convert egrid (m) -> e/D and evaluate f(Re,e/D)
    eD_vec = egrid./D;
    f = f_handle(Re0*ones(size(eD_vec)), eD_vec);

    ok = isfinite(f) & isreal(f) & (f>0);
    if nnz(ok) < 12
        viol_all(k) = 1;
        continue;
    end

    eok = egrid(ok);
    f_ok = f(ok);

    % deltaP at fixed Re0
    V0  = Re0 * mu_ref / (rho*D);
    DP  = f_ok .* ( rho .* (V0.^2) .* L ) ./ (2.*D);

    le  = log(eok);
    lDP = log(DP);
    if numel(lDP) >= smooth_w
        lDPs = movmean(lDP, smooth_w);
    else
        lDPs = lDP;
    end

    % s(e) = d log(ΔP) / d log(e)
    s = gradient(lDPs, le);

    % measure violation relative to band s >= -s_allow
    s_excess = min(s + s_allow, 0);  % <=0, how far below -s_allow
    bad      = (s_excess < 0);       % points with s < -s_allow

    frac = sum(bad) / max(1,numel(s));
    mag_raw = -mean(s_excess(bad));
    if isempty(mag_raw) || ~isfinite(mag_raw)
        mag_raw = 0;
    end
    mag = min(1, mag_raw / (s_allow + eps));

    viol_all(k) = min(1, 0.7*frac + 0.3*mag);
end

pen = mean(viol_all,'omitnan');
if ~isfinite(pen), pen = 1; end
end


function c3 = compute_C3_penalty_param(f_handle, phys)
% C3 in [0,1] (0 = best).
% Here: check deltaP vs mu at fixed V and e/D.

eD = 0.001;
mu = logspace(-4,-1,20).';
Re = (phys.rho*2.0*phys.D)./mu;

f  = f_handle(Re, eD*ones(size(Re)));
if any(~isfinite(f)) || any(~isreal(f)) || any(f<=0)
    c3 = 1;
    return;
end
DP = f .* (phys.rho*2.0^2*phys.L) ./ (2*phys.D);

dP = gradient(DP, mu);
% fraction of points where deltaP decreases more than small tolerance
tol = 0;  % we can make this >0 if we want to allow tiny non-monotonic bits
c3 = sum(dP < -tol) / numel(dP);
c3 = min(max(real(c3),0),1);
end


function c4 = compute_C4_penalty_param(f_handle, phys)
% C4 in [0,1] (0 = best).
% Check deltaP vs rho at fixed V and e/D, and γ ≈ d ln deltaP / d ln rho close to 1.

eD  = 0.001;
rho = logspace(2,4,20).';
Re  = (rho*2.0*phys.D)./phys.mu_ref;

f   = f_handle(Re, eD*ones(size(Re)));
if any(~isfinite(f)) || any(~isreal(f)) || any(f<=0)
    c4 = 1;
    return;
end

DP = f .* (rho*2.0^2*phys.L) ./ (2*phys.D);
ok = isfinite(DP) & (DP>0);
if nnz(ok) < 8
    c4 = 1;
    return;
end
rho_ok = rho(ok);
DP_ok  = DP(ok);

lf   = log(DP_ok);
lRho = log(rho_ok);
gamma = gradient(lf, lRho);   % ~1 ideally

c4 = mean(abs(gamma - 1.0));
c4 = min(max(real(c4),0),1);
end


function f_handle = assemble_model_handle_param(geneTemplatesX, p)
    numGenes = numel(geneTemplatesX);
    nTheta   = numGenes + 1;
    theta    = p(1:nTheta);
    C        = p(nTheta+1:end);
    C        = C(:).';
    f_handle = @(Re,eD) eval_model_param(geneTemplatesX, theta, C, [Re(:), eD(:)]);
end

function individ_updated = update_individual_constants_in_pop(individ, const_new)
%UPDATE_INDIVIDUAL_CONSTANTS_IN_POP
%   individ   : 1 x numGenes cell array of gene strings, e.g.
%               g(g(f(f(b(x1,[6.896276]),x2), ... )),x2)
%   const_new : vector of optimized constants Cbest, in the same global
%               order as const0 from encode_numeric_constants_all.
%
%   Rewrites each ERC [c_old] -> [c_new] in sequence.
%   If the number of ERCs and the length of const_new disagree, this
%   function throws an error.

    if ~iscell(individ)
        error('update_individual_constants_in_pop: individ must be a cell array of gene strings.');
    end

    const_new = real(const_new(:));   % ensure column, drop any tiny imag parts
    nNew      = numel(const_new);
    k         = 1;                    % global constant index

    individ_updated = individ;

    % Match ERC-style bracketed numeric literals: [3], [3.14], [-1.2e-3], etc.
    constPattern = '\[([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\]';

    for g = 1:numel(individ)
        s = individ{g};
        if ~(ischar(s) && ~isempty(s))
            continue;
        end

        out     = '';
        lastIdx = 1;

        [startIdx, endIdx] = regexp(s, constPattern, 'start','end');

        for m = 1:numel(startIdx)
            
            out = [out, s(lastIdx:startIdx(m)-1)];

            if k > nNew
                error(['update_individual_constants_in_pop: ', ...
                       'ran out of const_new entries while rewriting gene %d. ', ...
                       'Found more [c] literals than constants (used %d, have %d).'], ...
                       g, k-1, nNew);
            end

            newVal = const_new(k);
            out    = [out, '[', num2str(newVal,16), ']'];  %#ok<AGROW>
            k      = k + 1;

            lastIdx = endIdx(m) + 1;
        end

        
        out = [out, s(lastIdx:end)];
        individ_updated{g} = out;
        expr_g = individ_updated{g};
        
        n_open_par = sum(expr_g == '(');
        n_clos_par = sum(expr_g == ')');
        n_open_sq  = sum(expr_g == '[');
        n_clos_sq  = sum(expr_g == ']');

        if n_open_par ~= n_clos_par || n_open_sq ~= n_clos_sq
            error('update_individual_constants_in_pop: unbalanced brackets in gene %d: "%s"', ...
                g, expr_g);
        end

    end

    used = k - 1;
    if used ~= nNew
        error(['update_individual_constants_in_pop: mismatch between ERC count ', ...
               'and const_new length. Rewrote %d constants, but const_new has %d.'], ...
               used, nNew);
    end
end

