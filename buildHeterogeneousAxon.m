function [par, info] = buildHeterogeneousAxon(par0, lengths_um, varargin)
%BUILDHETEROGENEOUSAXON  Axon with per-internode (heterogeneous) lengths, clamped.
%
%   [par, info] = BUILDHETEROGENEOUSAXON(par0, lengths_um, ...)
%
%   Step-2 analogue of buildClampedAxon. Instead of a single scalar internode
%   length it installs a VECTOR of per-internode lengths (lengths_um, 1 x N, µm)
%   while holding the same strict-isolation clamp: axon radius and g-ratio are
%   reset to their reference values, the node geometry is unchanged, and the
%   paranodal periaxonal seal is re-applied (per internode) after the rebuild.
%   The number of internodes N is set by numel(lengths_um); the total internode
%   span is sum(lengths_um). Generate lengths_um with sampleInternodeLengths and
%   its 'TotalLength' option so the sum (and hence the clamp) is exactly matched
%   to the homogeneous comparator (see demo_step2_smoke / the step-2 plan).
%
%   At CV_L = 0 (all lengths equal) this reproduces buildClampedAxon exactly --
%   that is the builder's acceptance test.
%
%   Inputs:
%       par0        - baseline parameter struct (e.g. Carcamo2017CortexAxon()).
%       lengths_um  - 1 x N vector of per-internode TOTAL lengths (µm), > 0.
%
%   Name/value options:
%       'Nseg'                 - explicit internode segment count (overrides the
%                                resolution rule below).
%       'ResolveBy'            - 'max' (default) or 'mean': which length sets the
%                                global segment count. 'max' resolves the LONGEST
%                                internode at the step-1 segment length
%                                (L0/nseg0), guaranteeing every internode is at
%                                least that well resolved. See the plan §4/§9.
%       'MinSegments'          - floor on segment count (default 11).
%       'MaxSegments'          - cap on segment count (default [], none). Use a
%                                small cap ONLY for fast smoke tests.
%       'PreserveParanode'     - re-apply the paranodal seal (default true).
%       'ParanodeLength_um'    - paranode physical length held constant (1.9).
%       'ParanodeSealWidth_nm' - paranodal periaxonal width (0.012255632).
%
%   Output:
%       par   - parameter struct for the heterogeneous axon.
%       info  - struct: .nIntn, .nNode, .nIntSeg, .lengths_um (realized),
%               .meanL_um, .totalLength_um, .CV_L_actual,
%               .minSegLength_um, .maxSegLength_um.
%
%   See also: buildClampedAxon, sampleInternodeLengths, UpdateInternodeLength,
%             UpdateInternodePeriaxonalSpaceWidth, sweepHetCVEnergy.

nsegOpt    = getOption(varargin, 'Nseg', []);
resolveBy  = getOption(varargin, 'ResolveBy', 'max');
minSeg     = getOption(varargin, 'MinSegments', 11);
maxSeg     = getOption(varargin, 'MaxSegments', []);
preservePn = getOption(varargin, 'PreserveParanode', true);
pnLen_um   = getOption(varargin, 'ParanodeLength_um', 1.9);
sealW_nm   = getOption(varargin, 'ParanodeSealWidth_nm', 0.012255632);

% -------- validate --------
lengths_um = lengths_um(:).';                 % row vector
nIntn = numel(lengths_um);
if nIntn < 1
    error('buildHeterogeneousAxon:N', 'lengths_um must have at least one element.');
end
if any(~isfinite(lengths_um)) || any(lengths_um <= 0)
    error('buildHeterogeneousAxon:lengths', 'All internode lengths must be finite and > 0.');
end
nNode = nIntn + 1;

L0    = par0.intn.geo.length.value.ref;       % baseline internode length (µm)
nseg0 = par0.geo.nintseg;                      % baseline internode segment count

% -------- global internode segment count --------
% Hold the step-1 segment length (L0/nseg0) for the reference internode, so
% segment length = L_i/nseg <= L0/nseg0 for every internode when ResolveBy='max'.
if isempty(nsegOpt)
    switch lower(resolveBy)
        case 'max',  Lref = max(lengths_um);
        case 'mean', Lref = mean(lengths_um);
        otherwise,   error('buildHeterogeneousAxon:resolveBy', 'ResolveBy must be ''max'' or ''mean''.');
    end
    nseg = round(nseg0 * Lref / L0);
else
    nseg = round(nsegOpt);
end
nseg = max(minSeg, nseg);
if ~isempty(maxSeg), nseg = min(maxSeg, nseg); end

% -------- rebuild geometry (clamps radius & g-ratio to ref; wipes the seal) --------
par = UpdateNumberOfNodes(par0, nNode, 'reset', 'max', false);
par = UpdateNumberOfInternodeSegments(par, nseg);

% -------- install per-internode total lengths (case 1b: 1 x #internodes) --------
par = UpdateInternodeLength(par, lengths_um);

% -------- re-apply the paranodal seal, per internode (single combined call) --------
% nPn varies with internode length so the physical paranode length stays ~pnLen_um.
if preservePn
    allI = zeros(1, 0);
    allS = zeros(1, 0);
    for ii = 1:nIntn
        nPn  = max(1, round(pnLen_um * nseg / lengths_um(ii)));
        nPn  = min(nPn, floor(nseg / 2));
        segs = [1:nPn, (nseg - nPn + 1):nseg];
        allI = [allI, ii * ones(1, numel(segs))];   %#ok<AGROW>
        allS = [allS, segs];                          %#ok<AGROW>
    end
    par = UpdateInternodePeriaxonalSpaceWidth(par, sealW_nm * ones(1, numel(allS)), ...
                                              allI, allS, 'min');
end

% -------- info --------
if nargout >= 2
    Lvec = par.intn.geo.length.value.vec(:).';   % realized per-internode totals (µm)
    info = struct('nIntn', nIntn, 'nNode', nNode, 'nIntSeg', nseg, ...
                  'lengths_um', Lvec, ...
                  'meanL_um', mean(Lvec), ...
                  'totalLength_um', sum(Lvec), ...
                  'CV_L_actual', std(Lvec) / mean(Lvec), ...
                  'minSegLength_um', min(lengths_um) / nseg, ...
                  'maxSegLength_um', max(lengths_um) / nseg);
end
end


% ===================== local helper =====================
function val = getOption(args, name, default)
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name)
        val = args{k + 1};
    end
end
end
