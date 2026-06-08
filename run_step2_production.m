function OUT = run_step2_production(varargin)
%RUN_STEP2_PRODUCTION  Full Step-2 run: heterogeneity sweep, precision, trade-off, figures.
%
%   OUT = RUN_STEP2_PRODUCTION('Name', value, ...)
%
%   Produces the final step-2 result set and saves CSVs + figures to an output
%   folder. Three blocks (each toggleable):
%     (1) CV_L sweep at a representative mean -> speed, energy, precision vs CV_L.
%     (2) Trade-off figure (H2): each realization mapped to normalized
%         (speed, 1/jitter, 1/energy); heterogeneous cloud vs the homogeneous
%         (CV_L=0) point, with Pareto-nondominated realizations flagged.
%     (3) Homogeneous vs Poisson across mean internode length (figure-1 style),
%         with total length clamped (N varies with mean) -- speed & energy.
%
%   Geometry is matched between the deterministic and precision sweeps (same
%   BaseSeed + CVlevels), so realization r at a level is the SAME axon in both,
%   letting block (2) place each axon in the full 3-axis space.
%
%   Key Name/value options (all have sensible defaults):
%     'AxonFcn'        baseline axon (default @Carcamo2017CortexAxon)
%     'TotalLength_um' clamped total span (default baseline N*L0)
%     'MeanRep_um'     representative mean for blocks 1-2 (default baseline L0)
%     'CVlevels'       heterogeneity levels (default 0:0.1:1.0)
%     'MeanL_um'       means for block 3 (default [40 60 80 100 120])
%     'nReal'          deterministic realizations / level (default 40)
%     'nRealPrec'      precision realizations / level (default 15)
%     'nTrials'        noisy trials / realization (default 100)
%     'NoiseAmp_nA'    noise intensity -- CALIBRATE with test_step2 first (default 0.05)
%     'BaseSeed'       sampler/geometry seed base (default 7000)
%     'DoPrecision'    run the precision sweep (default true)
%     'DoAcrossMean'   run block 3 (default true)
%     'Parallel'       parfor for precision trials (default true)
%     'OutDir'         results folder (default <repo>/step2_results)
%     (builder options MaxSegments/ResolveBy/... are forwarded to the sweeps.)
%
%   Returns OUT with the sweep result structs and the output folder.
%   NOTE: at full settings the precision block is the heavy one; start with a
%   smaller nRealPrec/nTrials, and confirm resolution + noise with test_step2.

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

axonFcn  = getOpt(varargin, 'AxonFcn',        @Carcamo2017CortexAxon);
par0     = axonFcn();
L0       = par0.intn.geo.length.value.ref;
totalLen = getOpt(varargin, 'TotalLength_um', par0.geo.nintn * L0);
meanRep  = getOpt(varargin, 'MeanRep_um',     L0);
CVlevels = getOpt(varargin, 'CVlevels',       0:0.1:1.0);
meanLvec = getOpt(varargin, 'MeanL_um',       [40 60 80 100 120]);
nReal    = getOpt(varargin, 'nReal',          40);
nRealP   = getOpt(varargin, 'nRealPrec',      15);
nTrials  = getOpt(varargin, 'nTrials',        100);
noiseAmp = getOpt(varargin, 'NoiseAmp_nA',    0.05);
baseSeed = getOpt(varargin, 'BaseSeed',       7000);
doPrec   = getOpt(varargin, 'DoPrecision',    true);
doAcross = getOpt(varargin, 'DoAcrossMean',   true);
parOn    = getOpt(varargin, 'Parallel',       true);
outdir   = getOpt(varargin, 'OutDir',         fullfile(thisDir, 'step2_results'));
bopts    = builderOptions(varargin);
if ~exist(outdir, 'dir'), mkdir(outdir); end

sweepOpts = [{'AxonFcn', axonFcn, 'TotalLength_um', totalLen, 'BaseSeed', baseSeed}, bopts];

fprintf('\n===== STEP 2 PRODUCTION RUN =====\nOutput folder: %s\n', outdir);
fprintf('Clamp: total=%.0f um, meanRep=%.1f um, N=%d; CV_L in [%g..%g] (%d levels)\n', ...
        totalLen, meanRep, round(totalLen/meanRep), CVlevels(1), CVlevels(end), numel(CVlevels));

% =============================================================== Block 1
fprintf('\n[1] CV_L sweep at mean=%.1f um (deterministic: %d real)\n', meanRep, nReal);
Rc = sweepHetCVEnergy(meanRep, CVlevels, sweepOpts{:}, 'nReal', nReal, 'Plot', true, ...
        'SaveCsv', fullfile(outdir, 'cvenergy_vs_CVL.csv'));
saveFig(gcf, outdir, 'fig1_speed_energy_vs_CVL');

Rp = [];
if doPrec
    fprintf('[1b] precision sweep at mean=%.1f um (%d real x %d trials)\n', meanRep, nRealP, nTrials);
    Rp = sweepHetPrecision(meanRep, CVlevels, sweepOpts{:}, 'nReal', nRealP, 'nTrials', nTrials, ...
            'NoiseAmp_nA', noiseAmp, 'Parallel', parOn, ...
            'SaveCsv', fullfile(outdir, 'precision_vs_CVL.csv'));
    sp_std = rowStd(Rp.sigma_per_mm);
    f = figure('Name', 'Precision vs heterogeneity', 'Color', 'w');
    errorbar(Rp.CV_L, Rp.sigma_per_mm_mean, sp_std, '-o', 'LineWidth', 1.3);
    xlabel('Internode-length CV_L'); ylabel('Timing jitter \sigma_t per mm (ms/mm)');
    title('Precision vs heterogeneity'); grid on;
    saveFig(f, outdir, 'fig2_precision_vs_CVL');
end

% =============================================================== Block 2
fprintf('[2] Trade-off figure (speed-precision-energy)\n');
try
    plotTradeoff(Rc, Rp, outdir);
catch e
    warning('run_step2:tradeoff', 'Trade-off plot failed: %s', e.message);
end

% =============================================================== Block 3
acrossTable = [];
if doAcross
    fprintf('[3] Homogeneous vs Poisson across mean L (clamped total, deterministic)\n');
    nM  = numel(meanLvec);
    cvH = nan(nM,1); cvP = nan(nM,1); cvPs = nan(nM,1);
    eH  = nan(nM,1); eP  = nan(nM,1); ePs  = nan(nM,1);
    for k = 1:nM
        Rk = sweepHetCVEnergy(meanLvec(k), [0 1], sweepOpts{:}, 'nReal', nReal, ...
                'Plot', false, 'Verbose', false);
        cvH(k) = Rk.cv_mean(1);     cvP(k) = Rk.cv_mean(2);     cvPs(k) = Rk.cv_std(2);
        eH(k)  = Rk.charge_mean(1); eP(k)  = Rk.charge_mean(2); ePs(k)  = Rk.charge_std(2);
        fprintf('    mean=%5.1f um: CV  homog=%.3f  Poisson=%.3f | E homog=%.3f  Poisson=%.3f\n', ...
                meanLvec(k), cvH(k), cvP(k), eH(k), eP(k));
    end
    f = figure('Name', 'Homogeneous vs Poisson across mean L', 'Color', 'w');
    subplot(1,2,1);
        plot(meanLvec, cvH, '-o', 'LineWidth', 1.6); hold on;
        errorbar(meanLvec, cvP, cvPs, '--s', 'LineWidth', 1.6);
        xlabel('Mean internode length (\mum)'); ylabel('Conduction velocity (m/s)');
        title('Speed'); legend('Homogeneous (CV_L=0)', 'Poisson (CV_L=1)', 'Location', 'best'); grid on;
    subplot(1,2,2);
        plot(meanLvec, eH, '-o', 'LineWidth', 1.6); hold on;
        errorbar(meanLvec, eP, ePs, '--s', 'LineWidth', 1.6);
        xlabel('Mean internode length (\mum)'); ylabel('Na^+ charge / AP / mm (pC/mm)');
        title('Metabolic cost'); legend('Homogeneous', 'Poisson', 'Location', 'best'); grid on;
    saveFig(f, outdir, 'fig3_homog_vs_poisson_vs_meanL');
    acrossTable = table(meanLvec(:), cvH, cvP, cvPs, eH, eP, ePs, ...
        'VariableNames', {'meanL_um','cv_homog','cv_poisson','cv_poisson_std', ...
                          'E_homog','E_poisson','E_poisson_std'});
    try, writetable(acrossTable, fullfile(outdir, 'across_meanL_homog_vs_poisson.csv')); catch, end
end

OUT = struct('Rc', Rc, 'Rp', Rp, 'acrossTable', acrossTable, 'outdir', outdir, ...
             'config', struct('totalLen', totalLen, 'meanRep', meanRep, 'CVlevels', CVlevels, ...
                              'nReal', nReal, 'nRealPrec', nRealP, 'nTrials', nTrials, ...
                              'noiseAmp_nA', noiseAmp, 'baseSeed', baseSeed));
fprintf('\n===== DONE. Figures + CSVs in %s =====\n', outdir);
end


% ======================= trade-off plot =======================
function plotTradeoff(Rc, Rp, outdir)
% Map each (CV_L, realization) axon to normalized (speed, 1/jitter, 1/energy),
% all "higher = better"; plot the heterogeneous cloud vs the homogeneous point.
hasP = ~isempty(Rp);
nr = Rc.nReal; if hasP, nr = min(nr, Rp.nReal); end

levmat = repmat(Rc.CV_L, 1, nr);
sp = reshape(Rc.cv(:, 1:nr),               [], 1);
en = reshape(Rc.charge_pC_per_mm(:, 1:nr), [], 1);
lv = reshape(levmat,                       [], 1);
if hasP
    pr = reshape(1 ./ Rp.sigma_per_mm(:, 1:nr), [], 1);   % precision = 1/jitter
    ok = isfinite(sp) & isfinite(en) & isfinite(pr) & (en > 0);
else
    pr = [];
    ok = isfinite(sp) & isfinite(en) & (en > 0);
end
sp = sp(ok); en = en(ok); lv = lv(ok);
spN  = nz(sp);
enN  = nz(1 ./ en);                       % energy efficiency
if hasP
    pr = pr(ok);
    prN = nz(pr);
else
    prN = [];
end

isHom = (lv == min(lv));                  % CV_L = 0 baseline

if hasP
    P  = [spN, prN, enN];
    nd = paretoFront(P);

    f = figure('Name', 'Trade-off: speed-precision-energy', 'Color', 'w');
    scatter3(spN, prN, enN, 26, lv, 'filled'); hold on;
    scatter3(spN(isHom), prN(isHom), enN(isHom), 160, 'k', 'p', 'filled');   % homogeneous (star)
    scatter3(spN(nd), prN(nd), enN(nd), 70, 'r', 'o', 'LineWidth', 1.4);     % Pareto (ring)
    xlabel('speed (norm)'); ylabel('1/jitter (norm)'); zlabel('1/energy (norm)');
    title('Heterogeneous cloud vs homogeneous (\bigstar);  Pareto-nondominated ringed');
    cb = colorbar; cb.Label.String = 'CV_L'; grid on; view(135, 22);
    saveFig(f, outdir, 'fig4_tradeoff_3D');

    g = figure('Name', 'Trade-off projections', 'Color', 'w');
    proj(1, spN, enN, lv, isHom, 'speed (norm)', '1/energy (norm)');
    proj(2, spN, prN, lv, isHom, 'speed (norm)', '1/jitter (norm)');
    proj(3, prN, enN, lv, isHom, '1/jitter (norm)', '1/energy (norm)');
    saveFig(g, outdir, 'fig4b_tradeoff_projections');
else
    f = figure('Name', 'Trade-off: speed vs energy', 'Color', 'w');
    proj(0, spN, enN, lv, isHom, 'speed (norm)', '1/energy (norm)');
    title('Speed vs energy (run with DoPrecision for the full 3-axis trade-off)');
    saveFig(f, outdir, 'fig4_tradeoff_speed_energy');
end
end

function proj(spIdx, x, y, lv, isHom, xl, yl)
if spIdx > 0, subplot(1, 3, spIdx); end
scatter(x, y, 22, lv, 'filled'); hold on;
scatter(x(isHom), y(isHom), 140, 'k', 'p', 'filled');
xlabel(xl); ylabel(yl); grid on;
cb = colorbar; cb.Label.String = 'CV_L';
end

function y = nz(x)
% min-max normalize to [0,1]
mn = min(x); rg = max(x) - mn;
if rg <= 0, y = zeros(size(x)); else, y = (x - mn) / rg; end
end

function nd = paretoFront(P)
% nd(i)=true if row i is not dominated (higher is better on every column).
m = size(P, 1); nd = true(m, 1);
for i = 1:m
    for j = 1:m
        if all(P(j, :) >= P(i, :)) && any(P(j, :) > P(i, :))
            nd(i) = false; break
        end
    end
end
end


% ======================= small helpers =======================
function s = rowStd(X)
s = nan(size(X, 1), 1);
for i = 1:size(X, 1)
    v = X(i, ~isnan(X(i, :)));
    if numel(v) > 1, s(i) = std(v); end
end
end

function saveFig(fig, outdir, name)
try
    saveas(fig, fullfile(outdir, [name '.png']));
    savefig(fig, fullfile(outdir, [name '.fig']));
catch e
    warning('run_step2:saveFig', '%s', e.message);
end
end

function opts = builderOptions(args)
names = {'Nseg','ResolveBy','MinSegments','MaxSegments','PreserveParanode', ...
         'ParanodeLength_um','ParanodeSealWidth_nm'};
opts = {};
for k = 1:2:numel(args) - 1
    if any(strcmpi(args{k}, names)), opts(end+1:end+2) = {args{k}, args{k+1}}; end %#ok<AGROW>
end
end

function val = getOpt(args, name, default)
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name), val = args{k + 1}; end
end
end
