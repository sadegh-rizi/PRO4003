function e = energyConsumption(CURRENTS, INTERNODE_LENGTH, varargin)
%ENERGYCONSUMPTION  Metabolic-cost proxy: Na+ charge per action potential per mm.
%
%   e = ENERGYCONSUMPTION(CURRENTS, INTERNODE_LENGTH, ...)
%
%   Computes the time-integrated absolute Na+ current per propagating action
%   potential per unit distance (Na+ charge per AP per mm) over the interior
%   nodes of the axon. This is the energy axis used in the speed-precision-
%   energy trade-off: the Na+ that enters during each spike must be pumped
%   back out by the Na+/K+-ATPase (3 Na+ per ATP), so Na+ influx is a direct
%   proxy for the metabolic cost of signalling (Attwell & Laughlin, 2001).
%
%   CURRENTS must come from Model(par, [], verbose, recordCurrents) with
%   recordCurrents >= 1. The required integral of |I| over time is computed
%   inside Model and stored in CURRENTS.charge (nC), so this function needs no
%   time step and no voltage trace.
%
%   Inputs:
%       CURRENTS         - struct from Model() with fields .charge (NxC, nC),
%                          .channames (1xC cell) and .nodeColumns.
%       INTERNODE_LENGTH - (1 x (nNodes-1)) internode lengths (mm) from Model().
%
%   Name/value options:
%       'Ion'            - substring identifying the ion channels to sum
%                          (case-insensitive). Default 'Na' (matches
%                          'Fast Na+' and 'Persistent Na+').
%       'nAP'            - number of propagating APs to normalise by. Default 1.
%       'NodeSegPerNode' - node segments per node (par.geo.nnodeseg). Default 1.
%                          If >1, segment charges are summed within each node.
%       'DiscardNodes'   - nodes discarded at each end (default ~10% of nNodes).
%       'NodeRange'      - explicit [first last] interior node range (overrides DiscardNodes).
%       'NodePositions'  - explicit cumulative node positions (mm) (overrides INTERNODE_LENGTH).
%       'MinNodes'       - minimum interior nodes required (default 3).
%
%   Output struct e:
%       .chargePerAPperMM   - Na+ charge per AP per mm (nC/mm)   <-- primary readout.
%       .naChargePerNode    - per-node Na+ charge over the whole run (nC), interior nodes.
%       .totalNaCharge      - total interior Na+ charge over the whole run (nC).
%       .totalNaChargePerAP - total interior Na+ charge per AP (nC).
%       .interiorLength_mm  - axon length spanned by the interior nodes (mm).
%       .nAP                - number of APs used for normalisation.
%       .channelsUsed       - names of the channels summed as Na+.
%       .channelIdx         - their indices into CURRENTS.channames.
%       .naIonsPerAPperMM   - Na+ ions per AP per mm (charge / elementary charge).
%       .atpPerAPperMM      - ATP molecules per AP per mm (naIons / 3).
%       .units              - struct of unit strings.
%
%   See also: Model, conductionSpeed, measureAxon, interiorNodeIndices.

% -------- options --------
ion       = getOption(varargin, 'Ion', 'Na');
nAP       = getOption(varargin, 'nAP', 1);
nnodeseg  = getOption(varargin, 'NodeSegPerNode', 1);
discardN  = getOption(varargin, 'DiscardNodes', []);
nodeRng   = getOption(varargin, 'NodeRange', []);
nodePos   = getOption(varargin, 'NodePositions', []);
minNodes  = getOption(varargin, 'MinNodes', 3);

% -------- validate input --------
if ~isstruct(CURRENTS) || ~isfield(CURRENTS, 'charge') || isempty(CURRENTS.charge)
    error('energyConsumption:noCurrents', ...
        ['CURRENTS.charge is empty. Run Model with current recording on, e.g.\n', ...
         '    [V, IL, t, CURRENTS] = Model(par, [], true, 1);']);
end
if nAP <= 0
    error('energyConsumption:badnAP', 'nAP must be a positive number.');
end

% -------- identify Na+ channels by name --------
channames = CURRENTS.channames;
isNa = cellfun(@(s) ischar(s) && ~isempty(strfind(lower(s), lower(ion))), channames);
chIdx = find(isNa);
if isempty(chIdx)
    error('energyConsumption:noIon', ...
        'No channel name contains "%s". Available: %s', ion, strjoin(channames, ', '));
end

% -------- Na+ charge per node segment, then per node --------
naChargeSeg = sum(CURRENTS.charge(:, chIdx), 2);   % (nSeg x 1), nC
nnodeseg = round(nnodeseg);
if nnodeseg > 1
    if mod(numel(naChargeSeg), nnodeseg) ~= 0
        error('energyConsumption:badSeg', ...
            'Number of node segments (%d) is not divisible by NodeSegPerNode (%d).', ...
            numel(naChargeSeg), nnodeseg);
    end
    % Columns of CURRENTS.charge are ordered node-by-node with nnodeseg
    % consecutive segments per node, so reshape and sum down the segments.
    naChargePerNode = sum(reshape(naChargeSeg, nnodeseg, []), 1).';   % (nNodes x 1)
else
    naChargePerNode = naChargeSeg;
end
nNodes = numel(naChargePerNode);

% -------- node positions along the axon (mm) --------
if isempty(nodePos)
    if numel(INTERNODE_LENGTH) < nNodes - 1
        error('energyConsumption:badLengths', ...
            'INTERNODE_LENGTH must have at least nNodes-1 = %d elements.', nNodes - 1);
    end
    il = INTERNODE_LENGTH(:).';
    nodePos = [0, cumsum(il(1:nNodes-1))];
else
    nodePos = nodePos(:).';
end

% -------- interior nodes (same selection rule as conductionSpeed) --------
interior = interiorNodeIndices(nNodes, discardN, nodeRng, minNodes);

interiorCharge   = sum(naChargePerNode(interior));            % nC, all APs
interiorLength   = nodePos(interior(end)) - nodePos(interior(1));   % mm
chargePerAP      = interiorCharge / nAP;                      % nC per AP
if interiorLength > 0
    chargePerAPperMM = chargePerAP / interiorLength;          % nC / AP / mm
else
    chargePerAPperMM = NaN;
end

% -------- ion / ATP equivalents --------
eCharge          = 1.602176634e-19;          % C per elementary charge
chargeC          = chargePerAPperMM * 1e-9;  % nC -> C  (per AP per mm)
naIonsPerAPperMM = chargeC / eCharge;
atpPerAPperMM    = naIonsPerAPperMM / 3;      % 3 Na+ extruded per ATP

% -------- assemble output --------
e = struct();
e.chargePerAPperMM   = chargePerAPperMM;
e.naChargePerNode    = naChargePerNode(interior).';
e.totalNaCharge      = interiorCharge;
e.totalNaChargePerAP = chargePerAP;
e.interiorLength_mm  = interiorLength;
e.nAP                = nAP;
e.channelsUsed       = channames(chIdx);
e.channelIdx         = chIdx;
e.nodesUsed          = interior;
e.naIonsPerAPperMM   = naIonsPerAPperMM;
e.atpPerAPperMM      = atpPerAPperMM;
e.units = struct('charge', 'nC', 'chargePerAPperMM', 'nC/mm', ...
                 'naIonsPerAPperMM', 'ions/mm', 'atpPerAPperMM', 'ATP/mm');
end


% ===================== local helper =====================

function val = getOption(args, name, default)
%GETOPTION  Minimal case-insensitive name/value option reader.
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name)
        val = args{k + 1};
    end
end
end
