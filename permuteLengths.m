function Lp = permuteLengths(L, mode, seed)
%PERMUTELENGTHS  Reorder a fixed internode-length multiset (arrangement knob).
%
%   Lp = PERMUTELENGTHS(L, mode, seed)
%
%   Returns a permutation of L (same values, same mean/total/CV_L) in a chosen
%   spatial order, to isolate the *arrangement* effect from the *distribution*
%   effect (Talidou & Lefebvre 2025: conduction is path-dependent, with velocity
%   rising long->short and falling short->long).
%
%   mode:
%       'random'      - random permutation (default). Pass seed for reproducibility.
%       'ascending'   - sorted short->long (few, monotonic transitions).
%       'descending'  - sorted long->short.
%       'alternating' - interleave high/low to maximise long<->short transitions.
%       'clustered'   - all long internodes grouped together (= ascending here).
%
%   Use with buildHeterogeneousAxon to compare orders of the same multiset.

if nargin < 2 || isempty(mode), mode = 'random'; end

L = L(:).';
n = numel(L);

switch lower(mode)
    case 'random'
        if nargin >= 3 && ~isempty(seed)
            prev = rng; cleaner = onCleanup(@() rng(prev)); rng(seed); %#ok<NASGU>
        end
        Lp = L(randperm(n));

    case {'ascending', 'clustered'}
        Lp = sort(L, 'ascend');

    case 'descending'
        Lp = sort(L, 'descend');

    case 'alternating'
        s  = sort(L, 'descend');
        Lp = zeros(1, n);
        nHi = ceil(n / 2);
        Lp(1:2:n) = s(1:nHi);            % high values in odd slots
        Lp(2:2:n) = s(n:-1:nHi + 1);     % low values in even slots
        % adjacent slots alternate high/low -> many sharp transitions

    otherwise
        error('permuteLengths:mode', ...
              'mode must be random | ascending | descending | alternating | clustered.');
end
end
