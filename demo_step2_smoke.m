function R = demo_step2_smoke()
%DEMO_STEP2_SMOKE  Fast self-checking smoke test for the Step-2 heterogeneity code.
%
%   R = DEMO_STEP2_SMOKE()
%
%   Runs a tiny, fast version of the whole step-2 pipeline and ASSERTs the
%   invariants, so you can confirm the new code runs and is wired correctly
%   before launching production sweeps. It exercises:
%       [1] sampleInternodeLengths  - sum/mean/positivity invariants
%       [2] buildHeterogeneousAxon  - CV_L=0 reproduces buildClampedAxon (key test)
%       [3] one deterministic heterogeneous run (CV + energy, must propagate)
%       [4] sweepHetCVEnergy        - tiny 2-level x 2-real sweep
%       [5] sweepHetPrecision       - tiny 2-level x 2-real x 5-trial sweep
%       [6] permuteLengths          - arrangement preserves the multiset
%
%   It deliberately uses a SHORT axon (N=12) and a small segment cap
%   (MaxSegments) so it finishes quickly; those caps are for speed only and
%   must be removed for production accuracy (do a resolution convergence check).
%
%   Run from anywhere -- it adds the repo to the path itself.
%   On success it prints "ALL CHECKS PASSED" and returns the sweep results.

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

fprintf('\n=== STEP 2 SMOKE TEST ===\n');

% ---- tiny, fast configuration ----
par0   = Carcamo2017CortexAxon();
meanL  = par0.intn.geo.length.value.ref;     % 81.7 um baseline internode length
N      = 12;                                  % short axon (fast)
total  = N * meanL;                           % clamped total = N*meanL
maxSeg = 24;                                   % segment cap: SPEED ONLY (not production)
bopts  = {'MaxSegments', maxSeg};
tol    = 1e-6;

% ---------- [1] sampler invariants ----------
L = sampleInternodeLengths(N, meanL, 'CV', 0.4, 'TotalLength', total, 'Seed', 1);
assert(numel(L) == N,                         '[1] wrong length count');
assert(abs(sum(L)  - total) <= tol*total,     '[1] sum != clamped total');
assert(abs(mean(L) - meanL) <= tol*meanL,     '[1] mean != meanL');
assert(all(L > 0),                            '[1] non-positive internode length');
fprintf('[1] sampler OK: N=%d  sum=%.4f (target %.4f)  mean=%.4f  CV_L=%.3f\n', ...
        N, sum(L), total, mean(L), std(L)/mean(L));

% ---------- [2] CV_L = 0 must reproduce buildClampedAxon ----------
[parA, infoA] = buildClampedAxon(par0, meanL, total, bopts{:});
[parB, infoB] = buildHeterogeneousAxon(par0, meanL * ones(1, N), bopts{:});
assert(infoA.nIntn == infoB.nIntn, '[2] internode count mismatch');
assert(parA.geo.nintseg == parB.geo.nintseg, '[2] segment count mismatch');
mA = localMeasure(parA);
mB = localMeasure(parB);
relCV = abs(mA.cv - mB.cv) / max(abs(mA.cv), eps);
relE  = abs(mA.e  - mB.e ) / max(abs(mA.e),  eps);
fprintf('[2] CV_L=0 vs buildClampedAxon: dCV=%.2e  dE=%.2e  (expect ~0)\n', relCV, relE);
assert(mA.prop && mB.prop, '[2] baseline AP did not propagate');
assert(relCV < 1e-4 && relE < 1e-4, '[2] CV_L=0 does NOT reproduce buildClampedAxon');
fprintf('    builder equivalence OK\n');

% ---------- [3] one deterministic heterogeneous run ----------
parH = buildHeterogeneousAxon(par0, L, bopts{:});
mH = localMeasure(parH);
fprintf('[3] het run: CV=%.3f m/s  E=%.4g nC/AP/mm  propagated=%d\n', mH.cv, mH.e, mH.prop);
assert(mH.prop, '[3] AP did not propagate on heterogeneous axon');

% ---------- [4] tiny deterministic sweep ----------
Rc = sweepHetCVEnergy(meanL, [0 0.5], 'nReal', 2, 'TotalLength_um', total, ...
        'MaxSegments', maxSeg, 'Plot', false, 'Verbose', false);
assert(all(isfinite(Rc.cv_mean)), '[4] sweepHetCVEnergy produced non-finite CV');
fprintf('[4] sweepHetCVEnergy OK: CV(0)=%.3f  CV(0.5)=%.3f m/s | E(0)=%.3g  E(0.5)=%.3g pC/mm\n', ...
        Rc.cv_mean(1), Rc.cv_mean(2), Rc.charge_mean(1), Rc.charge_mean(2));

% ---------- [5] tiny precision sweep ----------
Rp = sweepHetPrecision(meanL, [0 0.5], 'nReal', 2, 'nTrials', 5, 'TotalLength_um', total, ...
        'MaxSegments', maxSeg, 'NoiseAmp_nA', 0.05, 'Parallel', false, 'Verbose', false);
assert(all(Rp.nValid(:) >= 1), '[5] some precision conditions had zero valid trials (raise Tmax_ms / lower noise)');
fprintf('[5] sweepHetPrecision OK: sigma_t(0)=%.4g ms  sigma_t(0.5)=%.4g ms  (det.latency std(0.5)=%.4g ms)\n', ...
        Rp.sigma_t_mean(1), Rp.sigma_t_mean(2), Rp.detLatency_std(2));

% ---------- [6] arrangement preserves the multiset ----------
for m = {'ascending','descending','alternating','random'}
    Lp = permuteLengths(L, m{1}, 11);
    assert(abs(sum(Lp) - sum(L)) < 1e-9, ['[6] permute changed total for ' m{1}]);
    assert(isequal(sort(Lp), sort(L)),    ['[6] permute changed multiset for ' m{1}]);
end
fprintf('[6] permuteLengths OK: all orders preserve the length multiset\n');

fprintf('=== ALL CHECKS PASSED ===\n\n');

R = struct('cvEnergy', Rc, 'precision', Rp, 'sampleLengths', L);
end


% ===================== local helper =====================
function m = localMeasure(par)
%LOCALMEASURE  One deterministic run -> CV (m/s) and Na+ charge/AP/mm (nC/mm).
[V, IL, t, C] = Model(par, [], false, 1);
dt = t(2) - t(1);
sp = conductionSpeed(V, IL, dt);
en = energyConsumption(C, IL, 'NodeSegPerNode', par.geo.nnodeseg);
m.cv   = sp.cv;
m.e    = en.chargePerAPperMM;
m.prop = sp.propagated && ~sp.failed;
end
