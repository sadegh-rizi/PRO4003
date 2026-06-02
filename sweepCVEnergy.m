function R = sweepCVEnergy(Lvalues_um, varargin)
%SWEEPCVENERGY  Deterministic CV and Na+ energy across internode lengths (clamped total).
%
%   R = SWEEPCVENERGY(Lvalues_um, ...)
%
%   For each internode length L it builds a clamped-total-length axon
%   (N = round(totalLength/L) internodes, radius and g-ratio held, paranodal
%   seal preserved; see buildClampedAxon), runs the DETERMINISTIC model once,
%   and measures over a fixed physical window:
%       * conduction velocity (m/s) via conductionSpeed (threshold crossing), and
%       * Na+ charge per AP per mm via energyConsumption,
%   reported as charge (pC/mm) and ATP (mol per AP per mm, Q/3F).
%
%   No trials are needed: CV and energy are deterministic.
%
%   Name/value options:
%       'AxonFcn'        - baseline axon builder (default @Carcamo2017CortexAxon).
%       'TotalLength_um' - clamped total internode span. Default = baseline N*L0.
%       'Vcross'         - CV threshold-crossing voltage, mV (default -20).
%       'Method'         - 'regression' (default) | 'twopoint'.
%       'WindowFrac'     - [lo hi] fraction of total length defining the fixed
%                          measurement window (default [0.2 0.8]).
%       'Verbose'        - print a line per L (default true).
%       'Plot'           - draw CV-vs-L and energy-vs-L (default true).
%       'SaveCsv'        - path to write a results CSV (default [], none).
%       'Tmax_ms'        - override sim duration (ms); raise it if distal nodes
%                          don't fire (propagated = false). Default keeps the axon's.
%       'Dt_us'          - override the time step (us); coarser = faster.
%       (builder options ScaleSegments/MinSegments/MaxSegments/PreserveParanode
%        are forwarded to buildClampedAxon.)
%
%   Output struct R (one row per L): .L_um (nominal), .L_actual_um, .nIntn,
%       .nNode, .totalLength_um, .cv_m_per_s, .r2, .propagated,
%       .charge_pC_per_mm, .atp_mol_per_mm, .windowLength_mm, plus .table.
%
%   See also: buildClampedAxon, conductionSpeed, energyConsumption, sweepPrecision.

axonFcn  = getOption(varargin, 'AxonFcn', @Carcamo2017CortexAxon);
totalLen = getOption(varargin, 'TotalLength_um', []);
vcross   = getOption(varargin, 'Vcross', -20);
method   = getOption(varargin, 'Method', 'regression');
winFrac  = getOption(varargin, 'WindowFrac', [0.2 0.8]);
verbose  = getOption(varargin, 'Verbose', true);
doPlot   = getOption(varargin, 'Plot', true);
saveCsv  = getOption(varargin, 'SaveCsv', []);
tmaxMs   = getOption(varargin, 'Tmax_ms', []);
dtUs     = getOption(varargin, 'Dt_us', []);
builderOpts = builderOptions(varargin);

FARADAY = 96485;                                   % C/mol

par0 = axonFcn();
if isempty(totalLen)
    totalLen = par0.geo.nintn * par0.intn.geo.length.value.ref;   % baseline total
end

Lvalues_um = Lvalues_um(:).';
nL = numel(Lvalues_um);

L_act = nan(nL,1); nIntn = nan(nL,1); nNode = nan(nL,1); totL = nan(nL,1);
cvArr = nan(nL,1); r2Arr = nan(nL,1); propArr = false(nL,1);
qpc   = nan(nL,1); atp = nan(nL,1); dwin = nan(nL,1);

if verbose
    fprintf('Clamped-total sweep: total = %.0f um, %d lengths\n', totalLen, nL);
    fprintf('%7s %6s %6s %10s %8s %14s\n','L[um]','nIntn','R^2','CV[m/s]','win[mm]','E[pC/AP/mm]');
end

for q = 1:nL
    L = Lvalues_um(q);
    try
        [par, info] = buildClampedAxon(par0, L, totalLen, builderOpts{:});
        if ~isempty(tmaxMs), par.sim.tmax.value = tmaxMs; par.sim.tmax.units = {1,'ms',1}; end
        if ~isempty(dtUs),   par.sim.dt.value   = dtUs;   par.sim.dt.units   = {1,'us',1};  end
        [V, IL, t, C] = Model(par, [], false, 1);     % deterministic, currents on
        dt  = t(2) - t(1);
        pos = [0, cumsum(IL(:).')];                    % mm, one per node
        Ltot_actual = pos(end);

        % Fixed physical measurement window -> node index range.
        x1 = winFrac(1) * Ltot_actual;
        x2 = winFrac(2) * Ltot_actual;
        f  = find(pos >= x1, 1, 'first');
        l  = find(pos <= x2, 1, 'last');
        if isempty(f) || isempty(l) || l - f < 2          % guard: widen if too narrow
            f = 2; l = numel(pos) - 1;
        end
        nodeRange = [f l];

        cv = conductionSpeed(V, IL, dt, 'Method', method, 'Vcross', vcross, ...
                             'NodeRange', nodeRange);
        en = energyConsumption(C, IL, 'NodeSegPerNode', par.geo.nnodeseg, ...
                               'NodeRange', nodeRange);

        L_act(q)   = info.L_um;   nIntn(q) = info.nIntn;   nNode(q) = info.nNode;
        totL(q)    = Ltot_actual;
        cvArr(q)   = cv.cv;       r2Arr(q) = cv.r2;        propArr(q) = cv.propagated && ~cv.failed;
        qpc(q)     = en.chargePerAPperMM * 1000;           % nC/mm -> pC/mm
        atp(q)     = (en.chargePerAPperMM * 1e-9) / (3 * FARADAY);   % mol ATP / AP / mm
        dwin(q)    = en.interiorLength_mm;
    catch err
        if verbose, fprintf('  L = %g um FAILED: %s\n', L, err.message); end
    end
    if verbose
        fprintf('%7.1f %6d %6.3f %10.3f %8.3f %14.4g\n', ...
                L, nIntn(q), r2Arr(q), cvArr(q), dwin(q), qpc(q));
    end
end

R = struct();
R.L_um = Lvalues_um(:);  R.L_actual_um = L_act;  R.nIntn = nIntn;  R.nNode = nNode;
R.totalLength_um = totL; R.cv_m_per_s = cvArr;   R.r2 = r2Arr;     R.propagated = propArr;
R.charge_pC_per_mm = qpc; R.atp_mol_per_mm = atp; R.windowLength_mm = dwin;
R.clampedTotal_um = totalLen;
try
    R.table = table(R.L_um, R.L_actual_um, R.nIntn, R.cv_m_per_s, R.r2, ...
                    R.charge_pC_per_mm, R.atp_mol_per_mm, R.totalLength_um, R.propagated, ...
        'VariableNames', {'L_um','L_actual_um','nIntn','cv_m_per_s','r2', ...
                          'charge_pC_per_mm','atp_mol_per_mm','totalLength_um','propagated'});
catch
    R.table = [];
end

if ~isempty(saveCsv) && ~isempty(R.table)
    try, writetable(R.table, saveCsv); if verbose, fprintf('Wrote %s\n', saveCsv); end
    catch e, warning('sweepCVEnergy:csv','%s', e.message); end
end

if doPlot
    try
        ok = ~isnan(cvArr);
        figure('Name','Clamped sweep: CV & energy','Color','w');
        subplot(1,2,1); plot(L_act(ok), cvArr(ok), '-o','LineWidth',1.3);
        xlabel('Internode length L (\mum)'); ylabel('Conduction velocity (m/s)');
        title('Speed vs L'); grid on;
        subplot(1,2,2); plot(L_act(ok), qpc(ok), '-o','LineWidth',1.3);
        xlabel('Internode length L (\mum)'); ylabel('Na^+ charge per AP per mm (pC/mm)');
        title('Energy vs L'); grid on;
    catch e, warning('sweepCVEnergy:plot','%s', e.message); end
end
end


% ===================== local helpers =====================

function opts = builderOptions(args)
% Pull only the buildClampedAxon-relevant name/value pairs out of varargin.
names = {'ScaleSegments','MinSegments','MaxSegments','PreserveParanode', ...
         'ParanodeLength_um','ParanodeSealWidth_nm'};
opts = {};
for k = 1:2:numel(args) - 1
    if any(strcmpi(args{k}, names))
        opts(end+1:end+2) = {args{k}, args{k+1}}; %#ok<AGROW>
    end
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
