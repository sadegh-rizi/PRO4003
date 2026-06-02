function cv = conductionSpeed(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, dt, varargin)
%CONDUCTIONSPEED  Conduction velocity of a propagating AP along a myelinated axon.
%
%   cv = CONDUCTIONSPEED(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, dt, ...)
%
%   Measures conduction velocity from the time at which the action potential
%   crosses a fixed voltage threshold at each node, using only interior nodes
%   (boundary internodes discarded). Two estimators are provided:
%
%     'regression' (default) - fit node distance (mm) vs. crossing time (ms)
%                              by least squares across all interior nodes;
%                              the slope is the conduction velocity. This is
%                              robust to per-node jitter and is the
%                              recommended readout for parameter sweeps.
%     'twopoint'             - distance / time between the first and last
%                              interior node, computed via the existing
%                              velocities() helper.
%
%   Because internode length is in mm and time in ms, the slope mm/ms equals
%   m/s directly (1 mm/ms = 1 m/s).
%
%   Inputs:
%       MEMBRANE_POTENTIAL - (T x N) node membrane potential (mV), as returned
%                            by Model() (one column per node site).
%       INTERNODE_LENGTH   - (1 x (N-1)) centre-to-centre internode lengths (mm),
%                            as returned by Model().
%       dt                 - sample interval (ms), e.g. TIME_VECTOR(2)-TIME_VECTOR(1).
%
%   Name/value options:
%       'Method'        - 'regression' (default) | 'twopoint'.
%       'Vcross'        - threshold-crossing voltage (mV). Default -20.
%       'DiscardNodes'  - nodes discarded at each end (default ~10% of N).
%       'NodeRange'     - explicit [first last] interior node range (overrides DiscardNodes).
%       'NodePositions' - explicit cumulative node positions (mm) (overrides INTERNODE_LENGTH).
%       'Interp'        - true (default) for sub-sample linear interpolation of the
%                         crossing time; false to use the nearest sample.
%       'MinNodes'      - minimum interior nodes that must fire (default 3).
%
%   Output struct cv:
%       .cv           - conduction velocity (m/s); NaN if the AP failed to propagate.
%       .method       - estimator used.
%       .vcross       - threshold used (mV).
%       .nodesUsed    - interior node indices that fired and were used.
%       .arrivalTimes - crossing time (ms) at each used node.
%       .positions    - cumulative distance (mm) of each used node.
%       .r2           - R^2 of the distance-vs-time fit ('regression' only; NaN otherwise).
%       .propagated   - true if every selected interior node crossed Vcross.
%       .failed       - true if fewer than MinNodes interior nodes fired.
%
%   See also: velocities, energyConsumption, measureAxon, interiorNodeIndices.

% -------- options --------
method   = getOption(varargin, 'Method', 'regression');
vcross   = getOption(varargin, 'Vcross', -20);
discardN = getOption(varargin, 'DiscardNodes', []);
nodeRng  = getOption(varargin, 'NodeRange', []);
nodePos  = getOption(varargin, 'NodePositions', []);
doInterp = getOption(varargin, 'Interp', true);
minNodes = getOption(varargin, 'MinNodes', 3);

N = size(MEMBRANE_POTENTIAL, 2);

% -------- node positions along the axon (mm) --------
if isempty(nodePos)
    if numel(INTERNODE_LENGTH) < N - 1
        error('conductionSpeed:badLengths', ...
            'INTERNODE_LENGTH must have at least N-1 = %d elements.', N - 1);
    end
    il = INTERNODE_LENGTH(:).';
    nodePos = [0, cumsum(il(1:N-1))];
else
    nodePos = nodePos(:).';
end

% -------- interior nodes --------
interior = interiorNodeIndices(N, discardN, nodeRng, minNodes);

% -------- crossing (arrival) time at each interior node --------
arrival = nan(1, numel(interior));
for q = 1:numel(interior)
    arrival(q) = arrivalTime(MEMBRANE_POTENTIAL(:, interior(q)), dt, vcross, doInterp);
end

fired      = ~isnan(arrival);
propagated = all(fired);

% Default / failure outputs.
cv = struct('cv', NaN, 'method', method, 'vcross', vcross, ...
            'nodesUsed', interior(fired), 'arrivalTimes', arrival(fired), ...
            'positions', nodePos(interior(fired)), 'r2', NaN, ...
            'propagated', propagated, 'failed', false);

if sum(fired) < minNodes
    cv.failed = true;
    return
end

t = arrival(fired);            % ms
x = nodePos(interior(fired));  % mm

switch lower(method)
    case 'twopoint'
        % First and last firing interior node, via the shipped velocities() helper.
        firstNode = interior(find(fired, 1, 'first'));
        lastNode  = interior(find(fired, 1, 'last'));
        cv.cv = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, dt, ...
                           [firstNode, lastNode], 'voltagecross', vcross);
        cv.r2 = NaN;

    case 'regression'
        % Distance (mm) = cv * time (ms) + b  ->  slope is the velocity (m/s).
        p = polyfit(t, x, 1);
        cv.cv = p(1);
        xfit  = polyval(p, t);
        ssres = sum((x - xfit).^2);
        sstot = sum((x - mean(x)).^2);
        if sstot > 0
            cv.r2 = 1 - ssres / sstot;
        else
            cv.r2 = NaN;
        end

    otherwise
        error('conductionSpeed:badMethod', ...
            'Unknown Method "%s" (use ''regression'' or ''twopoint'').', method);
end
end


% ===================== local helpers =====================
% (arrival-time logic lives in the shared AuxiliaryFunctions/arrivalTime.m)

function val = getOption(args, name, default)
%GETOPTION  Minimal case-insensitive name/value option reader.
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name)
        val = args{k + 1};
    end
end
end
