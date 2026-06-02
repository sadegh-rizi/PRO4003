function [par, info] = buildClampedAxon(par0, L_um, totalLength_um, varargin)
%BUILDCLAMPEDAXON  Axon at internode length L with the TOTAL length clamped.
%
%   [par, info] = BUILDCLAMPEDAXON(par0, L_um, totalLength_um, ...)
%
%   Implements the proposal's strict-isolation convention: the total axon
%   length is held (approximately) constant while the internode length L
%   varies, so the NUMBER of internodes changes as N = round(totalLength/L).
%   Axon radius and g-ratio are held at their reference values (the clamp),
%   and the paranodal periaxonal seal is re-applied after the rebuild.
%
%   Both sweepCVEnergy and sweepPrecision call this, so the two share exactly
%   the same geometry at every L.
%
%   Inputs:
%       par0           - baseline parameter struct (e.g. Carcamo2017CortexAxon()).
%                        Its reference values define radius, g-ratio, node length,
%                        channels, baseline internode length L0 and segment count.
%       L_um           - target internode length (micrometres).
%       totalLength_um - target total internode span (micrometres). N = round(total/L).
%
%   Name/value options:
%       'ScaleSegments'        - scale internode segment count with L so the
%                                segment LENGTH stays ~constant (default true).
%       'MinSegments'          - floor on internode segment count (default 11).
%       'MaxSegments'          - cap on segment count (default [], no cap).
%       'PreserveParanode'     - re-apply the paranodal seal (default true).
%       'ParanodeLength_um'    - paranode length held constant (default 1.9).
%       'ParanodeSealWidth_nm' - paranodal periaxonal width (default 0.012255632).
%
%   Output:
%       par   - parameter struct for the clamped axon.
%       info  - struct: .nIntn, .nNode, .nIntSeg, .L_um (actual), .totalLength_um (actual).
%
%   See also: sweepCVEnergy, sweepPrecision, UpdateNumberOfNodes.

scaleSeg  = getOption(varargin, 'ScaleSegments', true);
minSeg    = getOption(varargin, 'MinSegments', 11);
maxSeg    = getOption(varargin, 'MaxSegments', []);
preservePn = getOption(varargin, 'PreserveParanode', true);
pnLen_um   = getOption(varargin, 'ParanodeLength_um', 1.9);
sealW_nm   = getOption(varargin, 'ParanodeSealWidth_nm', 0.012255632);

L0    = par0.intn.geo.length.value.ref;   % baseline internode length (um)
nseg0 = par0.geo.nintseg;                  % baseline internode segment count

% Number of internodes needed to (approximately) hit the target total length.
nIntn = max(1, round(totalLength_um / L_um));
nNode = nIntn + 1;

% Internode segment count (resolution): hold segment length ~constant.
if scaleSeg
    nseg = round(nseg0 * L_um / L0);
else
    nseg = nseg0;
end
nseg = max(minSeg, nseg);
if ~isempty(maxSeg)
    nseg = min(maxSeg, nseg);
end

% Rebuild the axon at the new node count from reference values (this clamps
% radius and g-ratio to ref, and resets the periaxonal vector -> seal wiped).
% warnUser=false avoids the interactive 'reset?' prompt.
par = UpdateNumberOfNodes(par0, nNode, 'reset', 'max', false);

% Set the segment resolution, then the internode length.
par = UpdateNumberOfInternodeSegments(par, nseg);   % rebuilds seg arrays (peri uniform)
par = UpdateInternodeLength(par, L_um);             % total internode length = L_um

% Re-apply the paranodal seal that the rebuild wiped, holding the paranode
% physical length (~pnLen_um) constant across L.
if preservePn
    nPn  = max(1, round(pnLen_um * par.geo.nintseg / L_um));
    nPn  = min(nPn, floor(par.geo.nintseg / 2));
    segs = [1:nPn, (par.geo.nintseg - nPn + 1):par.geo.nintseg];
    par  = UpdateInternodePeriaxonalSpaceWidth(par, sealW_nm, [], segs, 'min');
end

if nargout >= 2
    info = struct('nIntn', nIntn, 'nNode', nNode, 'nIntSeg', par.geo.nintseg, ...
                  'L_um', par.intn.geo.length.value.ref, ...
                  'totalLength_um', nIntn * par.intn.geo.length.value.ref);
end
end


function val = getOption(args, name, default)
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name)
        val = args{k + 1};
    end
end
end
