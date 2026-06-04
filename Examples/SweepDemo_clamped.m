% SweepDemo_clamped.m
% Demo of the clamped-total-length sweeps: conduction velocity, Na+ energy,
% and timing precision against internode length on ONE figure.
%
% Total axon length is held (approximately) constant while L varies, so the
% number of internodes changes as N = round(total/L) - the proposal convention.
% sweepCVEnergy and sweepPrecision share buildClampedAxon, so the axes are
% measured on the same geometry family.
%
%   >> SweepDemo_clamped
%
% Speed notes (why these settings):
%   * CV + energy are deterministic (one run per L), so they get a FINE, LOW grid
%     that reaches down to ~20 um - that is where the conduction-velocity peak
%     sits, so a grid starting at 40 um would miss it.
%   * Precision is the slow step (nTrials runs per L). It is monotone-ish, so a
%     COARSE grid of a few lengths is enough.
%   * tmax = 3 ms: at this axon's CV (~2.5 m/s) the AP reaches the distal node in
%     ~1.4 ms, so 3 ms is plenty - 8 ms just simulated empty time.
% For publication-quality precision use nTrials >= 100 (and split L across runs).

thisDir  = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(thisDir)));        % put the repo on the path (this session)

% ---- demo settings ----
L_cv    = 20:5:120;      % fine, low grid for CV + energy (deterministic, cheap)
L_prec  = 40:20:120;     % coarse grid for precision (noisy, expensive)
tmax    = 3;             % ms - AP reaches the distal node by ~1.4 ms; 3 ms is plenty
dtus    = 0.25;          % us
nTrials = 50;            % demo only; use >= 100 for a stable precision estimate
sigma   = 0.001;         % nA noise intensity - TUNE so jitter is sub-millisecond
L0      = 81.7;          % cortical baseline internode length (um), reference line

% ---- deterministic: conduction velocity + energy (one run per L) ----
fprintf('=== CV + energy (deterministic, %d lengths) ===\n', numel(L_cv));
R1 = sweepCVEnergy(L_cv, 'Tmax_ms', tmax, 'Dt_us', dtus, 'Plot', false, ...
                   'SaveCsv', fullfile(thisDir, 'cv_energy_demo.csv'));

% ---- noisy: timing precision (nTrials per L) ----
fprintf('\n=== precision (%d lengths x %d trials - the slow part) ===\n', numel(L_prec), nTrials);
R2 = sweepPrecision(L_prec, 'Tmax_ms', tmax, 'Dt_us', dtus, 'nTrials', nTrials, 'BaseSeed',1000,...
                    'NoiseAmp_nA', sigma, 'SaveCsv', fullfile(thisDir, 'precision_demo_2.csv'));

% ---- sanity checks: did the AP reach the measurement points in time? ----
if any(~R1.propagated)
    warning('SweepDemo:cv', 'CV/energy: %d length(s) did not propagate - raise tmax.', sum(~R1.propagated));
end
if any(R2.nValid < nTrials)
    warning('SweepDemo:prec', 'precision: distal node missed on some trials - raise tmax.');
end

% ---- combined figure: speed / energy / precision vs L ----
figure('Name', 'Speed / energy / precision vs L (clamped total)', 'Color', 'w');

subplot(1, 3, 1);
plot(R1.L_actual_um, R1.cv_m_per_s, '-o', 'LineWidth', 1.3); hold on;
yl = ylim; plot([L0 L0], yl, 'k--'); ylim(yl);
xlabel('Internode length L (\mum)'); ylabel('Conduction velocity (m/s)');
title('Speed'); grid on;

subplot(1, 3, 2);
plot(R1.L_actual_um, R1.charge_pC_per_mm, '-o', 'LineWidth', 1.3); hold on;
yl = ylim; plot([L0 L0], yl, 'k--'); ylim(yl);
xlabel('Internode length L (\mum)'); ylabel('Na^+ charge per AP per mm (pC/mm)');
title('Energy'); grid on;

subplot(1, 3, 3);
plot(R2.L_actual_um, R2.sigma_per_mm, '-o', 'LineWidth', 1.3); hold on;
yl = ylim; plot([L0 L0], yl, 'k--'); ylim(yl);
xlabel('Internode length L (\mum)'); ylabel('Timing jitter \sigma_t / d (ms/mm)');
title('Precision'); grid on;

% ---- objective-specific optima (per the three curves) ----
[~, is] = max(R1.cv_m_per_s);              % fastest
[~, ie] = min(R1.charge_pC_per_mm);        % cheapest (lowest Na+ per mm)
[~, ip] = min(R2.sigma_per_mm);            % most precise (lowest jitter)
fprintf('\nObjective optima over the tested grid:\n');
fprintf('  L_speed     (max CV)      ~ %g um\n', R1.L_actual_um(is));
fprintf('  L_energy    (min cost)    ~ %g um\n', R1.L_actual_um(ie));
fprintf('  L_precision (min jitter)  ~ %g um\n', R2.L_actual_um(ip));
fprintf(['If energy and precision both bottom out at the long-L end while speed peaks\n', ...
         'in the middle, that is H1: no single L optimises all three.\n']);
