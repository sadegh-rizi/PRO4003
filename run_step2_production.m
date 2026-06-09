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
%     'DoCVEnergy'     run the deterministic speed/energy sweep (default true)
%     'DoPrecision'    run the precision sweep (default true)
%     'DoAcrossMean'   run block 3 (default true)
%     'Parallel'       parfor for precision trials (default true)
%     'OutDir'         results folder (default <repo>/step2_results)
%     'FastPass'       apply SAFE speed-ups -- nTrials=30, auto-trimmed Tmax...

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
noiseAmp = getOpt(varargin, 'NoiseAmp_nA',    0.02);
baseSeed = getOpt(varargin, 'BaseSeed',       7000);
doCV     = getOpt(varargin, 'DoCVEnergy',     true); % <-- NEW TOGGLE
doPrec   = getOpt(varargin, 'DoPrecision',    true);
doAcross = getOpt(varargin, 'DoAcrossMean',   true);
parOn    = getOpt(varargin, 'Parallel',       true);
outdir   = getOpt(varargin, 'OutDir',         fullfile(thisDir, 'step2_results'));
fastPass = getOpt(varargin, 'FastPass',       false);
tmaxMs   = getOpt(varargin, 'Tmax_ms',        []);
dtUs     = getOpt(varargin, 'Dt_us',          []);
bopts    = builderOptions(varargin);
if ~exist(outdir, 'dir'), mkdir(outdir); end

if fastPass
    if ~hasOpt(varargin, 'nTrials'), nTrials = 30; end
    if ~hasOpt(varargin, 'Tmax_ms')
        tmaxMs = autoTmax(par0, meanRep, totalLen, CVlevels, bopts, baseSeed);
    end
    poolN = 0;
    if parOn
        g = gcp('nocreate');
        if isempty(g)
            try, g = parpool; catch e, warning('run_step2:pool', 'Could not start a pool: %s', e.message); g = []; end
        end
        if ~isempty(g), poolN = g.NumWorkers; end
    end
end

sweepOpts = [{'AxonFcn', axonFcn, 'TotalLength_um', totalLen, 'BaseSeed', baseSeed}, bopts];
if ~isempty(tmaxMs), sweepOpts = [sweepOpts, {'Tmax_ms', tmaxMs}]; end
if ~isempty(dtUs),   sweepOpts = [sweepOpts, {'Dt_us',   dtUs}];   end

fprintf('\n===== STEP 2 PRODUCTION RUN =====\nOutput folder: %s\n', outdir);

% =============================================================== Block 1: Deterministic
Rc = [];
if doCV
    fprintf('\n[1] CV_L sweep at mean=%.1f um (deterministic: %d real)\n', meanRep, nReal);
    Rc = sweepHetCVEnergy(meanRep, CVlevels, sweepOpts{:}, 'nReal', nReal, 'Plot', true, ...
            'SaveCsv', fullfile(outdir, 'cvenergy_vs_CVL.csv'));
    saveFig(gcf, outdir, 'fig1_speed_energy_vs_CVL');
else
    fprintf('\n[1] Deterministic Speed/Energy sweep SKIPPED per user request.\n');
end

% =============================================================== Block 1b: Precision
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

% =============================================================== Block 2: Trade-off
if doPrec && doCV
    fprintf('[2] Trade-off figure (speed-precision-energy)\n');
    try, plotTradeoff(Rc, Rp, outdir); catch e, warning('run_step2:tradeoff', '%s', e.message); end
else
    fprintf('[2] 3D Trade-off plot skipped (requires both Speed and Precision sweeps to be active).\n');
end

% =============================================================== Block 3: Across Mean
acrossTable = [];
if doAcross
    fprintf('[3] Homogeneous vs Poisson across mean L (clamped total, deterministic)\n');
    nM  = numel(meanLvec);
    cvH = nan(nM,1); cvP = nan(nM,1); cvPs = nan(nM,1);
    eH  = nan(nM,1); eP  = nan(nM,1); ePs  = nan(nM,1);
    for k = 1:nM
        Rk = sweepHetCVEnergy(meanLvec(k), [0 1], sweepOpts{:}, 'nReal', nReal, 'Plot', false, 'Verbose', false);
        cvH(k) = Rk.cv_mean(1);     cvP(k) = Rk.cv_mean(2);     cvPs(k) = Rk.cv_std(2);
        eH(k)  = Rk.charge_mean(1); eP(k)  = Rk.charge_mean(2); ePs(k)  = Rk.charge_std(2);
    end
    f = figure('Name', 'Homogeneous vs Poisson across mean L', 'Color', 'w');
    subplot(1,2,1); plot(meanLvec, cvH, '-o', 'LineWidth', 1.6); hold on; errorbar(meanLvec, cvP, cvPs, '--s', 'LineWidth', 1.6); grid on;
    subplot(1,2,2); plot(meanLvec, eH, '-o', 'LineWidth', 1.6); hold on; errorbar(meanLvec, eP, ePs, '--s', 'LineWidth', 1.6); grid on;
    saveFig(f, outdir, 'fig3_homog_vs_poisson_vs_meanL');
    acrossTable = table(meanLvec(:), cvH, cvP, cvPs, eH, eP, ePs, 'VariableNames', {'meanL_um','cv_homog','cv_poisson','cv_poisson_std','E_homog','E_poisson','E_poisson_std'});
    try, writetable(acrossTable, fullfile(outdir, 'across_meanL_homog_vs_poisson.csv')); catch, end
end

OUT = struct('Rc', Rc, 'Rp', Rp, 'acrossTable', acrossTable, 'outdir', outdir);
ends