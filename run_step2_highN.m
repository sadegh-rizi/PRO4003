%RUN_STEP2_HIGHN  Re-run the Step-2 heterogeneity sweep with more realizations
%   and more noisy trials than step2_example.m used.
%
%   WHY: precision_vs_CVL.csv currently has nRealPrec = 10 and nTrials = 25
%   (see step2_example.m), while the deterministic CV/energy sweep used
%   nReal = 30-40. That mismatch is why the "beats homogeneous" count in
%   fig_reach_precision.png (4 of ~48) has wide, noisy confidence bands, and
%   it's well short of the 300 trials the Methods section of the manuscript
%   currently claims. This script reruns BOTH blocks together, with a single
%   unified CV_L grid and matched BaseSeed, so realization r at level CV_L is
%   the SAME axon geometry in the CV/energy sweep and the precision sweep --
%   which is what lets fig4 (3D trade-off) and fig5 (dominance/reach tests)
%   pair them correctly.
%
%   STEP 1 (do this first): run the PILOT block below (small grid, small
%   nReal/nTrials) and time it. Model runs don't scale linearly with nTrials
%   because of parfor overhead, so estimate from a real pilot, not by
%   multiplying the old runtime by the ratio of realizations.
%
%   STEP 2: once you know the per-run cost, size FULL.nReal / FULL.nTrials to
%   whatever your machine/time budget allows and run the FULL block. SaveCsv
%   checkpoints the table after every CV_L level, so a crash partway through
%   doesn't lose earlier levels -- just re-point OutDir/SaveCsv and resume from
%   the next level if that happens.

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

RUN_PILOT = true;   % <-- set false once you've timed the pilot and are ready for FULL
RUN_FULL  = false;  % <-- set true to run the full production sweep

% Unified CV_L grid: includes 0 (homogeneous baseline, sanity check) and 0.6,
% which nothing has tested under noise yet -- that's the gap between "100%
% propagated" at CV_L=0.6 (deterministic) and "0% valid" at CV_L=0.7 (noisy)
% flagged in the manuscript placeholder for the conduction-failure result.
CVlevels = 0:0.1:0.9;

% ===================================================================
% PILOT — small and fast, just to measure wall-clock time per (level x
% realization x trial) before committing to the full run.
% ===================================================================
if RUN_PILOT
    tic;
    Rpilot = run_step2_production( ...
        'AxonFcn',        @Carcamo2017CortexAxon, ...
        'TotalLength_um', 4085, ...
        'MeanRep_um',     81.7, ...
        'CVlevels',       [0 0.3 0.6 0.9], ...  % 4 levels only, spanning the range
        'MeanL_um',       [40 60 80 100 120], ...
        'nReal',          5, ...    % deterministic realizations/level (pilot)
        'nRealPrec',      5, ...    % precision realizations/level (pilot)
        'nTrials',        20, ...   % noisy trials/realization (pilot)
        'NoiseAmp_nA',    0.02, ... % keep whatever you already calibrated in test_step2 — don't re-guess this
        'BaseSeed',       7000, ...
        'DoPrecision',    true, ...
        'DoAcrossMean',   false, ...   % skip block 3 for the pilot, it's not what we're timing
        'Parallel',       true, ...
        'OutDir',         fullfile(thisDir, 'step2_results_pilot'));
    pilotSeconds = toc;
    fprintf('\nPILOT done in %.1f s for 4 levels x 5 real x 20 trials.\n', pilotSeconds);
    perRunSeconds = pilotSeconds / (4 * 5 * 20);
    fprintf('~%.3f s per (level x realization x trial) model run.\n', perRunSeconds);
    % Use this number to size FULL below before setting RUN_FULL = true.
    % Example: 10 levels x 40 real x 100 trials = 40000 runs.
    estFullSeconds = perRunSeconds * numel(CVlevels) * 40 * 100;
    fprintf('Full run at nReal=40, nTrials=100 would take roughly %.0f min (parfor scaling not accounted for).\n', ...
            estFullSeconds / 60);
end

% ===================================================================
% FULL — the actual re-run. Two nTrials options commented below:
%   100  matches sweepHetPrecision's own wrapper default (a reasonable step up
%        from 25 without exploding runtime).
%   300  matches what the manuscript's Methods section currently claims for
%        the homogeneous sweep -- use this if you want the heterogeneous
%        sweep to be run at the same trial count you're citing in the paper.
% ===================================================================
if RUN_FULL
    OUT = run_step2_production( ...
        'AxonFcn',        @Carcamo2017CortexAxon, ...
        'TotalLength_um', 4085, ...            % clamped total span = 50 x 81.7 um baseline
        'MeanRep_um',     81.7, ...            % cortical baseline mean internode length
        'CVlevels',       CVlevels, ...        % unified grid — same for CV/energy AND precision
        'MeanL_um',       [40 60 80 100 120], ...
        'nReal',          40, ...              % deterministic realizations/level (was 30-40, now consistent)
        'nRealPrec',      40, ...              % precision realizations/level (was 10 -- this is the main fix)
        'nTrials',        100, ...             % noisy trials/realization (was 25; see comment above re: 300)
        'NoiseAmp_nA',    0.02, ...            % CALIBRATED — see test_step2 Part C; don't change without re-checking
        'BaseSeed',       7000, ...            % matched between deterministic & precision sweeps -- keep this fixed
        'DoPrecision',    true, ...
        'DoAcrossMean',   true, ...
        'Parallel',       true, ...            % needs a parpool open, or it falls back to serial
        'OutDir',         fullfile(thisDir, 'step2_results_highN'));

    fprintf('\nFULL run complete. New CSVs/figures in %s\n', OUT.outdir);
    fprintf('Compare step2_results_highN/precision_vs_CVL.csv against the old\n');
    fprintf('step2_results/precision_vs_CVL.csv to see how much the CIs tightened\n');
    fprintf('and whether the "beats homogeneous" count in fig_reach_precision changes.\n');
end
