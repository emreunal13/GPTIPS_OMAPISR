function [symObj,cellsymObj] = gpmodel2sym(gp,ID,fastMode,useAlias)
%GPMODEL2SYM (Physics-Aware Version)
%   Creates a symbolic object for a model.
%   CRITICAL CHANGE: Uses stored weights (Nelder-Mead) instead of re-calculating SVD.

symObj = []; cellsymObj = [];

if ~gp.info.toolbox.symbolic
    error('The Symbolic Math Toolbox is required to use this function.');
end

if nargin < 2
    disp('Usage: SYMOBJ = GPMODEL2SYM(GP,ID)'); return;
end

if nargin < 3 || isempty(fastMode), fastMode = false; end
if nargin < 4 || isempty(useAlias), useAlias = false; end

% Handle ID selection ('best', 'valbest', etc)
if ischar(ID)
    if strcmpi(ID,'best')
        if isfield(gp.results,'best')
            ID = gp.results.best.index;
        else
            error('No "best" result field found.');
        end
    elseif strcmpi(ID,'valbest')
        if isfield(gp.results,'valbest')
            ID = gp.results.valbest.index;
        else
            error('No validation data or valbest found.');
        end
    elseif strcmpi(ID,'testbest')
        if isfield(gp.results,'testbest')
            ID = gp.results.testbest.index;
        else
            error('No test data or testbest found.');
        end
    end
end

if isnumeric(ID) && (ID < 1 || ID > gp.runcontrol.pop_size)
    error('ID must be between 1 and Population Size');
end

% Extract genes
if isstruct(ID) % User supplied struct
    model = ID;
    genes = model.genes;
    % For a struct, we often don't have stored theta easily unless passed.
    % Fallback to genes2gpmodel logic if needed, but usually ID is an index.
else
    % Get expression strings
    genes = tree2evalstr(gp.pop{ID},gp);
end


% RETRIEVE STORED WEIGHTS

usedStoredWeights = false;
theta = [];

% Check if this individual has stored weights from the run
if isnumeric(ID) && isfield(gp.fitness, 'returnvalues') && ...
   numel(gp.fitness.returnvalues) >= ID && ~isempty(gp.fitness.returnvalues{ID})
    
    theta = gp.fitness.returnvalues{ID};
    usedStoredWeights = true;
end

% Fallback: If no weights stored (e.g. fresh struct), re-calculate
if isempty(theta)
    warning('No stored weights found. Re-calculating using standard SVD (Physics constraints may be lost!).');
    % Evaluate genes on training data
    X = gp.userdata.xtrain;
    y = gp.userdata.ytrain;
    
    % Eval genes
    numGenes = numel(genes);
    G = ones(size(X,1), numGenes+1);
    
    % Safe Eval
    for i=1:numGenes
        str = genes{i};
        str = regexprep(str,'x(\d+)','X(:,$1)');
        G(:,i+1) = eval(str);
    end
    
    theta = pinv(G)*y;
end


% CONSTRUCT SYMBOLIC OBJECT
% Bias
term = sym(theta(1));
if useAlias
    digits(4);
    term = vpa(term);
end
symObj = term;
cellsymObj{1} = term;

% Genes
for i = 1:numel(genes)
    % Create symbolic gene
    geneStr = genes{i};
    
    % Create symbolic variable x1..xn
    
    w = theta(i+1);
    if useAlias, w=vpa(w); end
    
    % We construct the string "w * gene" and convert to sym
    % Caution: 'str2sym' or 'sym' usage varies by Matlab version.
    
    % Method: Construct full string then parse
    % This avoids variable definition issues
    
    % Add to total model
    sym_gene = str2sym(geneStr); 
    symObj = symObj + (w * sym_gene);
    
    cellsymObj{i+1} = (w * sym_gene);
end

% FAST MODE
if fastMode
    % symObj = simplify(symObj, 'Steps', 10); % Optional
end

end