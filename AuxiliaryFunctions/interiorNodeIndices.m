function idx = interiorNodeIndices(nNodes, discardNodes, nodeRange, minNodes)
%INTERIORNODEINDICES  Select interior nodes, discarding boundary nodes.
%   idx = INTERIORNODEINDICES(nNodes, discardNodes, nodeRange, minNodes)
%
%   Both conductionSpeed() and energyConsumption() use this helper so that the
%   speed and energy readouts are measured over the *same* interior nodes,
%   keeping the per-unit-distance quantities directly comparable.
%
%   Inputs:
%       nNodes       - total number of recorded node sites (columns of the
%                      membrane-potential matrix / rows of a per-node quantity).
%       discardNodes - number of nodes to discard at EACH end. If empty, a
%                      default of max(1, round(0.1*nNodes)) is used. Boundary
%                      internodes are discarded because the stimulus end and
%                      the sealed far end are not at conduction steady state.
%       nodeRange    - optional explicit [first last] node index range. If
%                      supplied (non-empty), it overrides discardNodes.
%       minNodes     - minimum number of interior nodes required (default 3).
%
%   Output:
%       idx          - row vector of interior node indices.
%
%   The discard count is clamped so that at least minNodes interior nodes are
%   always returned (falling back to all nodes for very short axons).

VariableDefault('discardNodes', []);
VariableDefault('nodeRange', []);
VariableDefault('minNodes', 3);

nNodes = round(nNodes);
if nNodes < 1
    error('interiorNodeIndices:badN', 'nNodes must be a positive integer.');
end
minNodes = max(1, min(minNodes, nNodes));

% Explicit range takes precedence.
if ~isempty(nodeRange)
    n1 = max(1, round(nodeRange(1)));
    n2 = min(nNodes, round(nodeRange(end)));
    if n2 < n1
        error('interiorNodeIndices:badRange', 'nodeRange must have first <= last.');
    end
    idx = n1:n2;
    return
end

% Default discard: ~10% of nodes at each end, at least one.
if isempty(discardNodes)
    discardNodes = max(1, round(0.1 * nNodes));
end
discardNodes = max(0, round(discardNodes));

% Clamp the discard so that at least minNodes interior nodes survive.
maxDiscard = floor((nNodes - minNodes) / 2);
if maxDiscard < 0
    maxDiscard = 0;
end
if discardNodes > maxDiscard
    discardNodes = maxDiscard;
end

idx = (1 + discardNodes):(nNodes - discardNodes);
if isempty(idx)
    idx = 1:nNodes;   % degenerate fallback (very short axon)
end
end
