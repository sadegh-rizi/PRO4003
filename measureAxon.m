function m = measureAxon(par, varargin)
%MEASUREAXON  Run the axon model once and return its speed and energy readouts.
%
%   m = MEASUREAXON(par, ...)
%
%   Convenience pipeline for parameter sweeps: it runs Model() a single time
%   with ionic-current recording enabled, then computes
%       * conduction velocity (m/s) with conductionSpeed(), and
%       * Na+ charge per AP per mm (nC/mm) with energyConsumption(),
%   over the same interior nodes so the two readouts are directly comparable.
%
%   Name/value options (forwarded to the underlying functions):
%       'Verbose'      - print Model progress (default false).
%       'Method'       - CV estimator, 'regression' (default) | 'twopoint'.
%       'Vcross'       - threshold-crossing voltage in mV (default -20).
%       'DiscardNodes' - nodes discarded at each end (default ~10% of nNodes).
%       'NodeRange'    - explicit [first last] interior node range.
%       'Ion'          - substring identifying Na+ channels (default 'Na').
%       'nAP'          - number of propagating APs (default 1).
%       'RecordLevel'  - 1 (default, integrated charge) or 2 (also full I trace).
%
%   Output struct m:
%       .cv_m_per_s              - conduction velocity (m/s).
%       .energy_nC_per_AP_per_mm - Na+ charge per AP per mm (nC/mm).
%       .atp_per_AP_per_mm       - ATP molecules per AP per mm.
%       .internodeLength_um      - homogeneous internode length, if available (um).
%       .propagated              - true if the AP propagated across the interior nodes.
%       .speed                   - full struct from conductionSpeed().
%       .energy                  - full struct from energyConsumption().
%       .dt_ms                   - simulation sample interval (ms).
%
%   See also: Model, conductionSpeed, energyConsumption.

verbose  = getOption(varargin, 'Verbose', false);
method   = getOption(varargin, 'Method', 'regression');
vcross   = getOption(varargin, 'Vcross', -20);
discardN = getOption(varargin, 'DiscardNodes', []);
nodeRng  = getOption(varargin, 'NodeRange', []);
ion      = getOption(varargin, 'Ion', 'Na');
nAP      = getOption(varargin, 'nAP', 1);
recLevel = getOption(varargin, 'RecordLevel', 1);

% Single model run with ionic-current recording.
[V, IL, t, CURRENTS] = Model(par, [], verbose, recLevel);
dt = t(2) - t(1);

% Speed and energy over the same interior window.
speed = conductionSpeed(V, IL, dt, 'Method', method, 'Vcross', vcross, ...
                        'DiscardNodes', discardN, 'NodeRange', nodeRng);
energy = energyConsumption(CURRENTS, IL, 'Ion', ion, 'nAP', nAP, ...
                           'NodeSegPerNode', par.geo.nnodeseg, ...
                           'DiscardNodes', discardN, 'NodeRange', nodeRng);

m = struct();
m.cv_m_per_s              = speed.cv;
m.energy_nC_per_AP_per_mm = energy.chargePerAPperMM;
m.atp_per_AP_per_mm       = energy.atpPerAPperMM;
if isfield(par, 'intn') && isfield(par.intn.geo.length.value, 'ref')
    m.internodeLength_um = par.intn.geo.length.value.ref;
else
    m.internodeLength_um = NaN;
end
m.propagated = speed.propagated && ~speed.failed;
m.speed      = speed;
m.energy     = energy;
m.dt_ms      = dt;
end


function val = getOption(args, name, default)
%GETOPTION  Minimal case-insensitive name/value option reader.
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name)
        val = args{k + 1};
    end
end
end
