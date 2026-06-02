# Running the sweeps — conduction velocity, energy, precision

A short guide for the team. All three axes are built from the same pieces, so
CV / energy and precision land on the same internode-length (L) axis and can be
compared directly.

## The building blocks

| Function | What it gives you |
|----------|-------------------|
| `buildClampedAxon(par0, L_um, total_um)` | a parameter struct `par` for an axon at internode length `L_um` with the **total length clamped** (N = round(total/L) internodes; radius and g-ratio held; paranodal seal preserved). |
| `Model(par, [], false, recordCurrents)` | `[V, IL, t, CURRENTS]` — node voltages `V` (mV, time × node), internode lengths `IL` (mm), time `t` (ms). Deterministic if `par` has no `.noise`; noisy if you set `par.noise`. `recordCurrents = 1` fills `CURRENTS` (Na⁺ charge); `0` skips it. |
| `conductionSpeed(V, IL, dt, ...)` | conduction velocity (`.cv` in m/s) from the −20 mV crossing. |
| `energyConsumption(CURRENTS, IL, ...)` | Na⁺ charge per AP per mm (`.chargePerAPperMM`, nC/mm) — convert to pC/mm and ATP below. |
| `arrivalTime(V(:,node), dt, vcross)` | spike arrival time (ms) at one node — the **same** crossing used by CV, so precision can't drift from it. |

Convention (decided): **total length is clamped**, so longer internodes mean
fewer nodes (N = total/L). Every readout is reported **per unit distance**.
`buildClampedAxon` does the geometry; you just pick `L`.

## Recipe 1 — CV + energy (deterministic, one run per L)

No trials needed: with noise off the model is deterministic.

Easiest — use the wrapper:

```matlab
R1 = sweepCVEnergy(30:10:120, 'Tmax_ms', 8, 'Dt_us', 0.25);
% R1.cv_m_per_s, R1.charge_pC_per_mm, R1.atp_mol_per_mm, R1.propagated, R1.table
```

What it does per length (manual version, for adapting):

```matlab
par0 = Carcamo2017CortexAxon();
Ltot = par0.geo.nintn * par0.intn.geo.length.value.ref;   % clamped total (~4085 um)

par = buildClampedAxon(par0, L, Ltot);
par.sim.tmax.value = 8;    par.sim.tmax.units = {1,'ms',1};
par.sim.dt.value   = 0.25; par.sim.dt.units   = {1,'us',1};

[V, IL, t, C] = Model(par, [], false, 1);          % deterministic, currents ON
dt = t(2) - t(1);

cv = conductionSpeed(V, IL, dt, 'Vcross', -20);    % cv.cv  (m/s)
en = energyConsumption(C, IL, 'NodeSegPerNode', par.geo.nnodeseg);

charge_pC_per_mm = en.chargePerAPperMM * 1000;                 % nC/mm -> pC/mm
atp_mol_per_mm   = (en.chargePerAPperMM * 1e-9) / (3 * 96485); % Q/(3F)
```

## Recipe 2 — precision (noisy, N trials per L)

Set `par.noise`, run many trials, take the spread of the distal arrival time.

Easiest — use the wrapper:

```matlab
R2 = sweepPrecision(30:10:120, 'NoiseAmp_nA', 0.002, 'nTrials', 100, ...
                    'Tmax_ms', 8, 'Dt_us', 0.25, 'SaveCsv', 'precision.csv');
% R2.sigma_t_ms, R2.sigma_per_mm, R2.nValid, R2.table
```

What it does per length (manual version — this is the precision module):

```matlab
par = buildClampedAxon(par0, L, Ltot);
par.sim.tmax.value = 8;    par.sim.tmax.units = {1,'ms',1};
par.sim.dt.value   = 0.25; par.sim.dt.units   = {1,'us',1};
par.noise.amp.value = 0.002;  par.noise.amp.units = {1,'nA',1};   % TUNE (see below)

arr = nan(1, nTrials);  distal = NaN;  d = NaN;  dt = NaN;
for k = 1:nTrials
    par.noise.seed = 1000 + k;                     % independent, reproducible
    [V, IL, t] = Model(par, [], false, 0);          % noisy run, currents OFF
    if k == 1
        dt  = t(2) - t(1);
        pos = [0, cumsum(IL(:).')];                 % node positions (mm)
        [~, distal] = min(abs(pos - 0.8*pos(end))); % distal node ~80% along
        d = pos(distal);                            % distance from stimulus (mm)
    end
    arr(k) = arrivalTime(V(:, distal), dt, vcross); % vcross = -20
end
good         = ~isnan(arr);
sigma_t      = std(arr(good));     % ms   (timing jitter)
sigma_per_mm = sigma_t / d;        % ms/mm
```

## Tuning and sanity checks (read this — it's where things go wrong)

- **Tune the noise amplitude.** After the dt-fix, `NoiseAmp_nA` is an intensity
  scaled by 1/√dt, not a per-step value. Too large and you get
  `Voltage has gone outside reasonable bounds` (the run diverged — its results
  are invalid). Start small (≈ 0.001–0.002 nA at dt = 0.25 µs), run one length,
  and raise it until `sigma_t` is sub-millisecond while the run stays stable.
- **Keep dt and noise amplitude fixed for the whole sweep.** The noise meaning
  is dt-dependent, so changing dt mid-study makes lengths non-comparable.
- **Check propagation.** CV/energy: if `propagated = false`, the AP didn't reach
  the far window — raise `Tmax_ms`. Precision: if `nValid < nTrials`, the distal
  node didn't fire on some trials — raise `Tmax_ms` (and/or lower the noise).
- **Stability vs speed.** dt = 0.25 µs is a good balance; dt = 1 µs is faster but
  can blow up on short internodes — if you see the bounds error with noise off,
  lower `Dt_us`.
- **Precision is the slow step** (nTrials × nL model runs). Run it in batches by
  calling `sweepPrecision` with subsets of the L grid; with `SaveCsv` it rewrites
  the CSV after every length, so partial results survive a crash.
- Expect a `Resetting the entire axon...` line per length — that's
  `buildClampedAxon` rebuilding the geometry, not an error.

## Putting the three together

Both sweeps call `buildClampedAxon`, so they share geometry at every L. For the
trade-off space: normalise each axis (z-score or min-max) and invert energy and
jitter so "higher = better" on all three, then locate `L_speed` (max CV),
`L_energy` (min charge/mm), `L_precision` (min σ/mm). A runnable end-to-end
example that plots all three vs L is in `Examples/SweepDemo_clamped.m`.

## Units cheat-sheet

| Axis | Quantity | Units |
|------|----------|-------|
| Speed | conduction velocity | m/s |
| Energy | Na⁺ charge per AP per mm | pC/mm (and ATP: mol per AP per mm) |
| Precision | distal-arrival jitter σ_t, and σ_t/d | ms, and ms/mm |
