function R = sweepInternodeLength(Lvalues_um, varargin)
%SWEEPINTERNODELENGTH  Sweep internode length and measure CV and Na+ energy.
%
%   R = SWEEPINTERNODELENGTH(Lvalues_um, ...)
%
%   Implements the Aim-1 deterministic sweep: vary homogeneous internode
%   length L while holding the NUMBER of internodes fixed, so total axon
%   length floats as N*L and the model stays structurally identical (same
%   node count, same Na+ sites) at every L. For each L it measures
%       * conduction velocity (m/s)            via conductionSpeed, and
%       * Na+ charge per AP per mm (nC/mm)      via energyConsumption,
%   over the same interior nodes (through measureAxon).
%
%   To avoid under-resolving long internodes, the number of internode
%   segments is scaled with L so the segment LENGTH stays ~constant at its
%   baseline value (a simple d_lambda-style resolution rule). Turn this off
%   with 'ScaleSegments', false to hold the segment count fixed instead.
%
%   Input:
%       Lvalues_um - vector of internode lengths to test (micrometres).
%
%   Name/value options:
%       'AxonFcn'       - handle building the baseline axon (default @Carcamo2017CortexAxon).
%       'Verbose'       - print a progress line per L (default true).
%       'Vcross'        - CV threshold-crossing voltage, mV (default -20).
%       'Method'        - CV estimator, 'regression' (default) | 'twopoint'.
%       'DiscardNodes'  - interior-node discard at each end (default ~10% of N).
%       'ScaleSegments' - scale segment count with L (default true).
%       'MinSegments'   - floor on internode segment count (default 11).
%       'MaxSegments'   - cap on internode segment count (default [], no cap).
%                         Capping speeds up large-L runs at the cost of resolution.
%       'Tmax_ms'       - override simulation duration (ms). Default [] keeps the
%                         axon's own tmax. Increase it if long axons fail to propagate.
%       'PreserveParanode' - restore the paranodal periaxonal seal after segments are
%                         rescaled (default true). The seal (a high-resistance periaxonal
%                         constriction at the internode ends) is essential to double-cable
%                         conduction; without it CV and the location of the CV peak change.
%       'ParanodeLength_um'    - paranode length held constant when restoring the seal (default 1.9).
%       'ParanodeSealWidth_nm' - paranodal periaxonal width (default 0.012255632, the Carcamo value).
%       'Plot'          - draw CV-vs-L and energy-vs-L (default true).
%       'SaveCsv'       - path to write a results CSV (default [], none).
%
%   Output struct R (one row per L; failed/non-propagating points are NaN):
%       .L_um, .cv_m_per_s, .energy_nC_per_AP_per_mm, .atp_per_AP_per_mm,
%       .propagated, .r2, .nIntSeg, .nIntn, .totalAxonLength_mm
%       .table   - the same as a MATLAB table (if available).
%       .baseline_L_um, .baseline_nIntSeg
%
%   Example:
%       R = sweepInternodeLength(20:10:200);
%       [~, k] = max(R.cv_m_per_s);  fprintf('CV peak near L = %g um\n', R.L_um(k));
%
%   See also: measureAxon, conductionSpeed, energyConsumption, UpdateInternodeLength.

axonFcn   = getOption(varargin, 'AxonFcn', @Carcamo2017CortexAxon);
verbose   = getOption(varargin, 'Verbose', true);
vcross    = getOption(varargin, 'Vcross', -20);
method    = getOption(varargin, 'Method', 'regression');
discardN  = getOption(varargin, 'DiscardNodes', []);
scaleSeg  = getOption(varargin, 'ScaleSegments', true);
minSeg    = getOption(varargin, 'MinSegments', 11);
maxSeg    = getOption(varargin, 'MaxSegments', []);
tmaxMs    = getOption(varargin, 'Tmax_ms', []);
doPlot    = getOption(varargin, 'Plot', true);
saveCsv   = getOption(varargin, 'SaveCsv', []);
preservePn = getOption(varargin, 'PreserveParanode', true);
pnLen_um   = getOption(varargin, 'ParanodeLength_um', 1.9);
sealW_nm   = getOption(varargin, 'ParanodeSealWidth_nm', 0.012255632);

Lvalues_um = Lvalues_um(:).';     % row
nL = numel(Lvalues_um);

% Baseline geometry (segment length we try to preserve, and fixed internode count).
base   = axonFcn();
L0     = base.intn.geo.length.value.ref;     % um
nseg0  = base.geo.nintseg;
nIntn  = base.geo.nintn;
nNode  = base.geo.nnode;
nodeL0 = base.node.geo.length.value.ref;     % um

% Preallocate.
cvArr   = nan(nL, 1);  enArr  = nan(nL, 1);  atpArr = nan(nL, 1);
propArr = false(nL, 1); r2Arr = nan(nL, 1);  segArr = nan(nL, 1);
totArr  = nan(nL, 1);

if verbose
    fprintf('Sweeping %d internode lengths (N = %d internodes held fixed)\n', nL, nIntn);
    fprintf('%8s %6s %10s %14s %10s %5s\n', 'L[um]', 'nseg', 'CV[m/s]', 'E[nC/AP/mm]', 'R^2', 'prop');
end

for q = 1:nL
    L = Lvalues_um(q);

    % --- segment count (resolution) ---
    if scaleSeg
        nseg = round(nseg0 * L / L0);
    else
        nseg = nseg0;
    end
    nseg = max(minSeg, nseg);
    if ~isempty(maxSeg)
        nseg = min(maxSeg, nseg);
    end

    % --- build this axon: fixed internode COUNT, new length & resolution ---
    try
        par = axonFcn();                 % baseline (carries the paranodal periaxonal seal)
        reapplySeal = false;
        if nseg ~= par.geo.nintseg
            % Re-segmenting rebuilds the periaxonal vector as a uniform value,
            % which DESTROYS the paranodal seal - so restore it afterwards.
            par = UpdateNumberOfInternodeSegments(par, nseg);
            reapplySeal = preservePn;
        end
        par = UpdateInternodeLength(par, L);               % total internode length = L um, spread over nseg
        if reapplySeal
            % Restore the high-resistance paranodal seal at both internode ends,
            % holding the paranode physical length (~pnLen_um) constant with L.
            nPn  = max(1, round(pnLen_um * par.geo.nintseg / L));
            nPn  = min(nPn, floor(par.geo.nintseg / 2));
            segs = [1:nPn, (par.geo.nintseg - nPn + 1):par.geo.nintseg];
            par  = UpdateInternodePeriaxonalSpaceWidth(par, sealW_nm, [], segs, 'min');
        end
        if ~isempty(tmaxMs)
            par.sim.tmax.value = tmaxMs;
            par.sim.tmax.units = {1, 'ms', 1};
        end

        m = measureAxon(par, 'Verbose', false, 'Vcross', vcross, ...
                        'Method', method, 'DiscardNodes', discardN);

        cvArr(q)   = m.cv_m_per_s;
        enArr(q)   = m.energy_nC_per_AP_per_mm;
        atpArr(q)  = m.atp_per_AP_per_mm;
        propArr(q) = m.propagated;
        r2Arr(q)   = m.speed.r2;
        segArr(q)  = nseg;
        totArr(q)  = (nIntn * L + nNode * nodeL0) / 1000;   % mm
    catch err
        if verbose
            fprintf('  L = %g um FAILED: %s\n', L, err.message);
        end
        segArr(q) = nseg;
        totArr(q) = (nIntn * L + nNode * nodeL0) / 1000;
    end

    if verbose
        fprintf('%8.1f %6d %10.3f %14.4g %10.4f %5d\n', ...
                L, segArr(q), cvArr(q), enArr(q), r2Arr(q), propArr(q));
    end
end

% --- assemble output ---
R = struct();
R.L_um                    = Lvalues_um(:);
R.cv_m_per_s              = cvArr;
R.energy_nC_per_AP_per_mm = enArr;
R.atp_per_AP_per_mm       = atpArr;
R.propagated              = propArr;
R.r2                      = r2Arr;
R.nIntSeg                 = segArr;
R.nIntn                   = nIntn;
R.totalAxonLength_mm      = totArr;
R.baseline_L_um           = L0;
R.baseline_nIntSeg        = nseg0;

% Optional MATLAB table.
try
    R.table = table(R.L_um, R.cv_m_per_s, R.energy_nC_per_AP_per_mm, ...
                    R.atp_per_AP_per_mm, R.r2, R.propagated, R.nIntSeg, R.totalAxonLength_mm, ...
        'VariableNames', {'L_um','cv_m_per_s','energy_nC_per_AP_per_mm', ...
                          'atp_per_AP_per_mm','r2','propagated','nIntSeg','totalAxonLength_mm'});
catch
    R.table = [];
end

% Optional CSV.
if ~isempty(saveCsv)
    writeCsv(saveCsv, R);
    if verbose, fprintf('Wrote results to %s\n', saveCsv); end
end

% Optional plot.
if doPlot
    try
        plotSweep(R, L0);
    catch err
        warning('sweepInternodeLength:plot', 'Could not plot: %s', err.message);
    end
end
end


% ===================== local helpers =====================

function plotSweep(R, L0)
ok = ~isnan(R.cv_m_per_s);
figure('Name', 'Internode-length sweep', 'Color', 'w');

subplot(1, 2, 1);
plot(R.L_um(ok), R.cv_m_per_s(ok), '-o', 'LineWidth', 1.3); hold on;
yl = ylim; plot([L0 L0], yl, 'k--'); ylim(yl);
xlabel('Internode length L (\mum)'); ylabel('Conduction velocity (m/s)');
title('Speed vs L'); grid on;
text(L0, yl(1) + 0.05*diff(yl), sprintf(' baseline %.0f', L0), 'Color', 'k');

subplot(1, 2, 2);
plot(R.L_um(ok), R.energy_nC_per_AP_per_mm(ok), '-o', 'LineWidth', 1.3); hold on;
yl = ylim; plot([L0 L0], yl, 'k--'); ylim(yl);
xlabel('Internode length L (\mum)'); ylabel('Na^+ charge per AP per mm (nC/mm)');
title('Energy vs L'); grid on;
end


function writeCsv(fname, R)
try
    if ~isempty(R.table)
        writetable(R.table, fname);
        return
    end
catch
end
% Manual fallback.
fid = fopen(fname, 'w');
if fid < 0, error('sweepInternodeLength:csv', 'Cannot open %s', fname); end
fprintf(fid, 'L_um,cv_m_per_s,energy_nC_per_AP_per_mm,atp_per_AP_per_mm,r2,propagated,nIntSeg,totalAxonLength_mm\n');
for i = 1:numel(R.L_um)
    fprintf(fid, '%g,%g,%g,%g,%g,%d,%g,%g\n', R.L_um(i), R.cv_m_per_s(i), ...
            R.energy_nC_per_AP_per_mm(i), R.atp_per_AP_per_mm(i), R.r2(i), ...
            R.propagated(i), R.nIntSeg(i), R.totalAxonLength_mm(i));
end
fclose(fid);
end


function val = getOption(args, name, default)
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name)
        val = args{k + 1};
    end
end
end
