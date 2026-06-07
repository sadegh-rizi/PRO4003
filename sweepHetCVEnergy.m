function R = sweepHetCVEnergy(meanL_um, CVlevels, varargin)
%SWEEPHETCVENERGY  Deterministic CV & Na+ energy vs heterogeneity level (fixed mean, clamped total).
%
%   R = SWEEPHETCVENERGY(meanL_um, CVlevels, ...)
%
%   Step-2 analogue of sweepCVEnergy. The mean internode length is held fixed;
%   the internode-length SPREAD (CV_L) is swept. For each CV_L level it draws
%   nReal heterogeneous axons (sampleInternodeLengths, sum clamped to N*meanL),
%   builds each with buildHeterogeneousAxon, runs the DETERMINISTIC model once,
%   and measures conduction velocity and Na+ charge/AP/mm over a fixed physical
%   window -- exactly the readouts of step 1, now as a distribution per level.
%   CV_L = 0 reproduces the homogeneous clamped axon (built-in baseline).
%
%   N is fixed at round(TotalLength/meanL), so every axon (and the homogeneous
%   baseline) shares N, total length, and mean -- only the distribution varies.
%
%   Name/value options:
%       'nReal'          - heterogeneous realizations per CV_L level (default 20).
%       'TotalLength_um' - clamped total internode span (default = baseline N*L0).
%       'AxonFcn'        - baseline axon builder (default @Carcamo2017CortexAxon).
%       'Order'          - arrangement of each draw: 'random' (default) |
%                          ascending | descending | alternating (see permuteLengths).
%       'SamplerMin'/'SamplerMax' - length bounds passed to the sampler (default 0/Inf).
%       'BaseSeed'       - sampler seed = BaseSeed + 1000*levelIdx + realization (default 7000).
%       'Vcross'         - CV threshold-crossing voltage, mV (default -20).
%       'Method'         - 'regression' (default) | 'twopoint'.
%       'WindowFrac'     - [lo hi] fraction of total length for the measurement window ([0.2 0.8]).
%       'Tmax_ms'/'Dt_us'- sim duration / time step overrides (speed / propagation).
%       'Verbose'        - print a line per realization (default true).
%       'Plot'           - errorbar CV & energy vs CV_L (default true).
%       'SaveCsv'        - path to write the flat results table (default []).
%       (builder options Nseg/ResolveBy/MinSegments/MaxSegments/PreserveParanode/...
%        are forwarded to buildHeterogeneousAxon.)
%
%   Output struct R:
%       .CV_L (levels), .nReal, .N, .meanL_um, .totalLength_um,
%       .cv (nLevels x nReal, m/s), .charge_pC_per_mm, .atp_mol_per_mm,
%       .realizedCV_L, .propagated (logical), and per-level summaries
%       .cv_mean/.cv_std, .charge_mean/.charge_std, .propagated_frac, plus .table.
%
%   See also: sweepCVEnergy, buildHeterogeneousAxon, sampleInternodeLengths, sweepHetPrecision.

axonFcn  = getOption(varargin, 'AxonFcn', @Carcamo2017CortexAxon);
nReal    = getOption(varargin, 'nReal', 20);
totalLen = getOption(varargin, 'TotalLength_um', []);
order    = getOption(varargin, 'Order', 'random');
sMin     = getOption(varargin, 'SamplerMin', 0);
sMax     = getOption(varargin, 'SamplerMax', inf);
baseSeed = getOption(varargin, 'BaseSeed', 7000);
vcross   = getOption(varargin, 'Vcross', -20);
method   = getOption(varargin, 'Method', 'regression');
winFrac  = getOption(varargin, 'WindowFrac', [0.2 0.8]);
tmaxMs   = getOption(varargin, 'Tmax_ms', []);
dtUs     = getOption(varargin, 'Dt_us', []);
verbose  = getOption(varargin, 'Verbose', true);
doPlot   = getOption(varargin, 'Plot', true);
saveCsv  = getOption(varargin, 'SaveCsv', []);
builderOpts = builderOptions(varargin);

FARADAY = 96485;                                  % C/mol

par0 = axonFcn();
if isempty(totalLen)
    totalLen = par0.geo.nintn * par0.intn.geo.length.value.ref;
end
N           = max(1, round(totalLen / meanL_um));
targetTotal = N * meanL_um;                       % match the homogeneous comparator exactly

CVlevels = CVlevels(:).';
nLev = numel(CVlevels);

cv   = nan(nLev, nReal);  qpc = nan(nLev, nReal);  atp = nan(nLev, nReal);
rCVL = nan(nLev, nReal);  prop = false(nLev, nReal);

if verbose
    fprintf('Het CV/energy sweep: meanL=%.1f um, N=%d, total=%.0f um, %d levels x %d real\n', ...
            meanL_um, N, targetTotal, nLev, nReal);
    fprintf('%6s %5s %9s %10s %8s %14s %5s\n', 'CV_L', 'real', 'CV_L(act)', 'CV[m/s]', 'R^2', 'E[pC/AP/mm]', 'prop');
end

for ci = 1:nLev
    CVL = CVlevels(ci);
    for r = 1:nReal
        try
            seed = baseSeed + 1000 * ci + r;
            L = sampleInternodeLengths(N, meanL_um, 'CV', CVL, 'TotalLength', targetTotal, ...
                                       'Min', sMin, 'Max', sMax, 'Seed', seed);
            if ~strcmpi(order, 'random'), L = permuteLengths(L, order); end

            par = buildHeterogeneousAxon(par0, L, builderOpts{:});
            if ~isempty(tmaxMs), par.sim.tmax.value = tmaxMs; par.sim.tmax.units = {1,'ms',1}; end
            if ~isempty(dtUs),   par.sim.dt.value   = dtUs;   par.sim.dt.units   = {1,'us',1};  end

            [V, IL, t, C] = Model(par, [], false, 1);            % deterministic, currents on
            dt  = t(2) - t(1);
            pos = [0, cumsum(IL(:).')];
            Ltot = pos(end);

            x1 = winFrac(1) * Ltot;  x2 = winFrac(2) * Ltot;
            f  = find(pos >= x1, 1, 'first');  l = find(pos <= x2, 1, 'last');
            if isempty(f) || isempty(l) || l - f < 2, f = 2; l = numel(pos) - 1; end
            nodeRange = [f l];

            spd = conductionSpeed(V, IL, dt, 'Method', method, 'Vcross', vcross, 'NodeRange', nodeRange);
            en  = energyConsumption(C, IL, 'NodeSegPerNode', par.geo.nnodeseg, 'NodeRange', nodeRange);

            cv(ci, r)   = spd.cv;
            qpc(ci, r)  = en.chargePerAPperMM * 1000;                 % nC/mm -> pC/mm
            atp(ci, r)  = (en.chargePerAPperMM * 1e-9) / (3 * FARADAY);
            rCVL(ci, r) = std(L) / mean(L);
            prop(ci, r) = spd.propagated && ~spd.failed;

            if verbose
                fprintf('%6.2f %5d %9.3f %10.3f %8.3f %14.4g %5d\n', ...
                        CVL, r, rCVL(ci, r), cv(ci, r), spd.r2, qpc(ci, r), prop(ci, r));
            end
        catch err
            if verbose, fprintf('  CV_L=%.2f real=%d FAILED: %s\n', CVL, r, err.message); end
        end
    end
end

% -------- assemble + per-level summaries --------
R = struct();
R.CV_L = CVlevels(:);  R.nReal = nReal;  R.N = N;
R.meanL_um = meanL_um;  R.totalLength_um = targetTotal;
R.cv = cv;  R.charge_pC_per_mm = qpc;  R.atp_mol_per_mm = atp;
R.realizedCV_L = rCVL;  R.propagated = prop;
R.cv_mean        = nanmeanRow(cv);     R.cv_std        = nanstdRow(cv);
R.charge_mean    = nanmeanRow(qpc);    R.charge_std    = nanstdRow(qpc);
R.propagated_frac = mean(prop, 2);

% flat table (one row per realization)
[CIg, Rg] = ndgrid(1:nLev, 1:nReal);
try
    R.table = table(CVlevels(CIg(:)).', Rg(:), rCVL(:), cv(:), qpc(:), atp(:), prop(:), ...
        'VariableNames', {'CV_L', 'realization', 'realizedCV_L', 'cv_m_per_s', ...
                          'charge_pC_per_mm', 'atp_mol_per_mm', 'propagated'});
catch
    R.table = [];
end
if ~isempty(saveCsv) && ~isempty(R.table)
    try, writetable(R.table, saveCsv); if verbose, fprintf('Wrote %s\n', saveCsv); end
    catch e, warning('sweepHetCVEnergy:csv', '%s', e.message); end
end

if doPlot
    try
        figure('Name', 'Het sweep: CV & energy vs CV_L', 'Color', 'w');
        subplot(1,2,1); errorbar(R.CV_L, R.cv_mean, R.cv_std, '-o', 'LineWidth', 1.3);
        xlabel('Internode-length CV_L'); ylabel('Conduction velocity (m/s)');
        title('Speed vs heterogeneity'); grid on;
        subplot(1,2,2); errorbar(R.CV_L, R.charge_mean, R.charge_std, '-o', 'LineWidth', 1.3);
        xlabel('Internode-length CV_L'); ylabel('Na^+ charge per AP per mm (pC/mm)');
        title('Energy vs heterogeneity'); grid on;
    catch e, warning('sweepHetCVEnergy:plot', '%s', e.message); end
end
end


% ===================== local helpers =====================
function m = nanmeanRow(X)
m = nan(size(X,1),1);
for i = 1:size(X,1), v = X(i,~isnan(X(i,:))); if ~isempty(v), m(i) = mean(v); end, end
end

function s = nanstdRow(X)
s = nan(size(X,1),1);
for i = 1:size(X,1), v = X(i,~isnan(X(i,:))); if numel(v) > 1, s(i) = std(v); end, end
end

function opts = builderOptions(args)
names = {'Nseg','ResolveBy','MinSegments','MaxSegments','PreserveParanode', ...
         'ParanodeLength_um','ParanodeSealWidth_nm'};
opts = {};
for k = 1:2:numel(args) - 1
    if any(strcmpi(args{k}, names)), opts(end+1:end+2) = {args{k}, args{k+1}}; end %#ok<AGROW>
end
end

function val = getOption(args, name, default)
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name), val = args{k + 1}; end
end
end
