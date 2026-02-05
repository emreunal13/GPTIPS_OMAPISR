function f = modelSym_to_fhandle(modelStr)
%MODELSYM_TO_FHANDLE Convert a modelSym string to a safe, vectorized function
%   f = modelSym_to_fhandle(modelStr)
%
%   Input:
%     modelStr - string of the model using variables x1 (Re) and x2 (e/D).
%                e.g. "1.5e-8*tanh(x1) + (6.5e-7/x1^(1/2) + 0.26*x2^(3/2) - 1.5e-8)/x2 + 0.011"
%
%   Output:
%     f - function handle f(Re,eD) that returns an array same size as inputs.
%         Invalid entries produce NaN.
%
% Notes:
%  - The routine attempts to make operators elementwise and to guard against
%    obvious numerical issues (log of nonpositive, division by zero, huge exponents).

% Basic sanity
if nargin<1 || isempty(modelStr)
    error('Provide a model string, e.g. modelSym output.');
end
if ~ischar(modelStr) && ~isstring(modelStr)
    error('modelStr must be a char or string.');
end
s = char(modelStr);

% 1) replace variable names x1->Re, x2->eD (word boundaries)
s = regexprep(s, '\<x1\>', 'Re');
s = regexprep(s, '\<x2\>', 'eD');

% 2) make operators elementwise:
%   - '^' -> '.^'
%   - '*' -> '.*' (but avoid already '.*')
%   - '/' -> './' (but avoid already './')
% Use lookbehind to avoid double-dot replacements.

% caret
s = regexprep(s, '(?<!\.)\^', '.^');

% multiplication: replace '*' that are NOT already '.*'
s = regexprep(s, '(?<!\.)\*', '.*');

% division: replace '/' that are NOT already './'
s = regexprep(s, '(?<!\.)/', './');

% For safety: convert integer power notation like x.^2 which might become .^. - already ok

% 3) ensure common functions are correctly named (tanh, exp, log, sqrt, sin, cos are fine)
% nothing to change here normally.

% 4) Build anonymous function string and convert using str2func
% We'll create a raw anonymous handle: @(Re,eD) ( <expr> )
anonExpr = ['@(Re,eD) ( ' s ' )'];

% Attempt to create function handle
try
    f_raw = str2func(anonExpr);
catch ME
    error('Failed to create function handle from model string. Expr:\n%s\nError: %s', anonExpr, ME.message);
end

% 5) Wrap in a safe evaluator that:
%    - broadcasts scalars to match arrays
%    - captures runtime errors and returns NaNs at those positions
%    - coerces nonreal/Inf to NaN
f = @(Re,eD) safe_eval(f_raw, Re, eD);

end

%% local helper
function out = safe_eval(f_raw, Re, eD)
    % ensure double arrays
    Re = double(Re); eD = double(eD);

    % broadcasting: allow scalar vs vector combinations
    if isscalar(Re) && ~isscalar(eD)
        Re = Re .* ones(size(eD));
    end
    if isscalar(eD) && ~isscalar(Re)
        eD = eD .* ones(size(Re));
    end

    % output container
    out = nan(size(Re));

    % mask where basic inputs are valid (Re>0 and eD>=0 usually)
    % allow eD==0 (some models may use it), but protect against negatives
    ok = isfinite(Re) & (Re>0) & isfinite(eD) & (eD>=0);

    if ~any(ok)
        return
    end

    % evaluate function in try/catch to catch domain errors
    try
        % evaluate only for valid positions to reduce runtime errors
        Re_sub = Re(ok);
        eD_sub = eD(ok);
        y = f_raw(Re_sub, eD_sub);  % f_raw should be vectorized
    catch
        % on error, return NaN for those positions
        out(ok) = NaN;
        return
    end

    % post-process y: convert Inf or complex to NaN, leave real finite values
    % if y is scalar but multiple ok positions, expand accordingly
    if isscalar(y) && numel(Re_sub)>1
        y = repmat(y, size(Re_sub));
    end

    % Make sure y shape matches
    if ~isequal(size(y), size(Re_sub))
        % try to reshape if possible
        try
            y = reshape(y, size(Re_sub));
        catch
            out(ok) = NaN;
            return
        end
    end

    % Guard values: non-real, non-finite -> NaN
    bad = ~isfinite(y) | ~isreal(y);
    y(bad) = NaN;

    % Put back into out
    out(ok) = y;
end
