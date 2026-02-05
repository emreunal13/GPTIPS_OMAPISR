function gp = gprandom(gp, seed)
%GPRANDOM  Set or initialize the random number generator for GPTIPS.
%   gp = gprandom(gp, seed)
%   If seed is omitted or empty, it uses the current clock time.
%
%   Example:
%     gp = gprandom(gp, 42);      % deterministic (recommended)
%     gp = gprandom(gp, 'shuffle'); % new random seed each run

% (c) Dominic Searson
%   GPTIPS 2 

% handle seed argument 
if nargin < 2 || isempty(seed)
    seed = sum(100 * clock); % default: time-dependent
elseif ischar(seed) && strcmpi(seed, 'shuffle')
    seed = sum(100 * clock);
elseif ~isnumeric(seed)
    error('Seed must be numeric or the string "shuffle".');
end

gp.info.PRNGseed = seed;

% create and assign stream
try
    s = RandStream('mt19937ar', 'Seed', gp.info.PRNGseed);
    RandStream.setGlobalStream(s);
catch
    % Fallback for older MATLAB
    rand('twister', gp.info.PRNGseed);
end

fprintf('[gprandom] RNG seed set to %d\n', gp.info.PRNGseed);

end
