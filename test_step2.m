function test_step2(runHeavyDiagnostics)
%TEST_STEP2  Validation + pre-production diagnostics for the Step-2 code.
%
%   test_step2()           runs everything (default).
%   test_step2(false)      runs only the fast asserted correctness checks (Part A).
%
%   Part A  (fast, ASSERTED): sampler invariants; CV_L=0 reproduces
%           buildClampedAxon; a heterogeneous axon propagates; the two sweeps run.
%   Part B  (heavier): resolution-convergence check at high CV_L -- compares CV
%           and energy at the default segment count vs 2x finer. This is the
%           check that tells you whether the energy-vs-CV_L rise you saw is real
%           or a under-resolution artifact (it should change <~1% if converged).
%   Part C  (heavier): noise-amplitude calibration -- prints distal-arrival jitter
%           for a few NoiseAmp_nA values so you can pick one giving sub-ms jitter
%           before the precision production run.
%
%   Run from anywhere; it adds the repo to the path itself.

if nargin < 1, runHeavyDiagnostics = true; end

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

fprintf('\n===== STEP 2 TEST & PRE-PRODUCTION DIAGNOSTICS =====\n');

par0  = Carcamo2017CortexAxon();
meanL = par0.intn.geo.length.value.ref;          % 81.7 um baseline

% ---------------------------------------------------------------- Part A
fprintf('\n[A] Correctness (fast, asserted)\n');
N = 12;  total = N * meanL;  maxSeg = 24;  bopts = {'MaxSegments', maxSeg};  tol = 1e-6;

L = sampleInternodeLengths(N, meanL, 'CV', 0.4, 'TotalLength', total, 'Seed', 1);
assert(numel(L) == N,                       'A1 sampler: wrong count');
assert(abs(sum(L)  - total) <= tol*total,   'A1 sampler: sum != clamped total');
assert(abs(mean(L) - meanL) <= tol*meanL,   'A1 sampler: mean != meanL');
assert(all(L > 0),                          'A1 sampler: non-positive length');
fprintf('    A1 sampler invariants OK (sum=%.3f, mean=%.3f, CV_L=%.3f)\n', sum(L), mean(L), std(L)/mean(L));

[pA, iA] = buildClampedAxon(par0, meanL, total, bopts{:});
[pB, iB] = buildHeterogeneousAxon(par0, meanL * ones(1, N), bopts{:});
assert(iA.nIntn == iB.nIntn && pA.geo.nintseg == pB.geo.nintseg, 'A2 geometry mismatch');
mA = measureOne(pA);  mB = measureOne(pB);
assert(mA.prop && mB.prop, 'A2 baseline AP did not propagate');
assert(abs(mA.cv - mB.cv)/abs(mA.cv) < 1e-4 && abs(mA.e - mB.e)/abs(mA.e) < 1e-4, ...
       'A2 CV_L=0 does NOT reproduce buildClampedAxon');
fprintf('    A2 builder equivalence OK (dCV=%.1e, dE=%.1e)\n', ...
        abs(mA.cv-mB.cv)/abs(mA.cv), abs(mA.e-mB.e)/abs(mA.e));

mH = measureOne(buildHeterogeneousAxon(par0, L, bopts{:}));
assert(mH.prop, 'A3 heterogeneous AP did not propagate');
fprintf('    A3 heterogeneous run OK (CV=%.3f m/s, E=%.4g nC/AP/mm)\n', mH.cv, mH.e);

Rc = sweepHetCVEnergy(meanL, [0 0.5], 'nReal', 2, 'TotalLength_um', total, ...
        'MaxSegments', maxSeg, 'Plot', false, 'Verbose', false);
Rp = sweepHetPrecision(meanL, [0 0.5], 'nReal', 2, 'nTrials', 5, 'TotalLength_um', total, ...
        'MaxSegments', maxSeg, 'NoiseAmp_nA', 0.05, 'Parallel', false, 'Verbose', false);
assert(all(isfinite(Rc.cv_mean)),  'A4 sweepHetCVEnergy produced non-finite CV');
assert(all(Rp.nValid(:) >= 1),     'A4 sweepHetPrecision had a zero-valid condition');
fprintf('    A4 sweeps run OK\n');
fprintf('[A] PASS\n');

if ~runHeavyDiagnostics
    fprintf('\n(Heavy diagnostics skipped.)\n===== TEST COMPLETE =====\n');
    return
end

% ---------------------------------------------------------------- Part B
fprintf('\n[B] Resolution convergence at high heterogeneity (CV_L=0.8)\n');
Nb = 20;  totalB = Nb * meanL;
Lh = sampleInternodeLengths(Nb, meanL, 'CV', 0.8, 'TotalLength', totalB, 'Seed', 101);
t0 = tic;
[pd, id] = buildHeterogeneousAxon(par0, Lh);                       % default ResolveBy='max'
md = measureOne(pd);
pf = buildHeterogeneousAxon(par0, Lh, 'Nseg', 2 * id.nIntSeg);     % 2x finer
mf = measureOne(pf);
dCV = abs(md.cv - mf.cv) / max(abs(mf.cv), eps) * 100;
dE  = abs(md.e  - mf.e ) / max(abs(mf.e ), eps) * 100;
fprintf('    nseg %d -> %d (%.1fs):  dCV=%.2f%%   dEnergy=%.2f%%\n', ...
        id.nIntSeg, 2*id.nIntSeg, toc(t0), dCV, dE);
if dE > 2 || dCV > 2
    fprintf('    *** NOT converged at default resolution: the energy/CV at high CV_L is\n');
    fprintf('        partly a discretisation artifact. Raise resolution (lower MaxSegments cap\n');
    fprintf('        / use ResolveBy=''max'' uncapped) for production. ***\n');
else
    fprintf('    Converged: default resolution is adequate at this CV_L.\n');
end

% ---------------------------------------------------------------- Part C
fprintf('\n[C] Noise-amplitude calibration (homogeneous, mean=%.1f um) -- aim for sub-ms jitter\n', meanL);
totalC = par0.geo.nintn * meanL;
amps = [0.02 0.05 0.1 0.2];
for a = amps
    Rn = sweepHetPrecision(meanL, 0, 'nReal', 1, 'nTrials', 15, 'TotalLength_um', totalC, ...
            'NoiseAmp_nA', a, 'Parallel', false, 'Verbose', false);
    fprintf('    amp=%.3f nA  ->  sigma_t=%.4g ms   (sigma/mm=%.4g ms/mm)\n', ...
            a, Rn.sigma_t_mean(1), Rn.sigma_per_mm_mean(1));
end
fprintf('    Pick the amp whose sigma_t is a small fraction of a ms and use it for production.\n');

fprintf('\n===== TEST COMPLETE =====\n');
end


% ===================== local helper =====================
function m = measureOne(par)
[V, IL, t, C] = Model(par, [], false, 1);
dt = t(2) - t(1);
sp = conductionSpeed(V, IL, dt);
en = energyConsumption(C, IL, 'NodeSegPerNode', par.geo.nnodeseg);
m.cv   = sp.cv;
m.e    = en.chargePerAPperMM;
m.prop = sp.propagated && ~sp.failed;
end
