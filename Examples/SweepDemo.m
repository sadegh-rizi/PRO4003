% SweepDemo - CV and Na+ energy cost across a sweep of internode lengths.
%
% Aim-1 deterministic sweep: vary homogeneous internode length L with the
% number of internodes held fixed (total axon length floats as N*L), and
% measure conduction velocity and Na+ charge per AP per mm at each L.
%
% Run from anywhere:
%   >> SweepDemo

% --- put the repository (and its subfolders) on the path for this session ---
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
addpath(genpath(repoRoot));

% --- choose the lengths to test (micrometres) ---
% Brackets the cortical baseline (~82 um). Start moderate: run time grows with
% L because the segment count is scaled up to keep spatial resolution constant.
Lvalues = [20 30 40 50 60 70 82 100 120 150 180 200];

% --- run the sweep (CV + energy at every L) ---
R = sweepInternodeLength(Lvalues, ...
        'AxonFcn',  @Carcamo2017CortexAxon, ...
        'Vcross',   -20, ...
        'Plot',     true, ...
        'SaveCsv',  fullfile(thisDir, 'SweepResults.csv'));

% --- locate the speed optimum and the cortical baseline ---
ok = ~isnan(R.cv_m_per_s);
[cvPeak, kPeak] = max(R.cv_m_per_s);
fprintf('\nConduction-velocity peak: %.2f m/s near L = %g um\n', cvPeak, R.L_um(kPeak));
fprintf('Cortical baseline internode length: %g um\n', R.baseline_L_um);

% Energy at the smallest vs largest L that propagated.
Lp = R.L_um(ok & R.propagated);
Ep = R.energy_nC_per_AP_per_mm(ok & R.propagated);
if numel(Lp) >= 2
    fprintf('Na+ energy/mm: %.3g nC/mm at L=%g um  ->  %.3g nC/mm at L=%g um\n', ...
            Ep(1), Lp(1), Ep(end), Lp(end));
end

% --- show the results table (if available) ---
if ~isempty(R.table)
    disp(' ');
    disp(R.table);
end

fprintf(['\nTip: any row with propagated = 0 means the AP did not reach all\n', ...
         'interior nodes within tmax - re-run those with a larger ''Tmax_ms''.\n']);
