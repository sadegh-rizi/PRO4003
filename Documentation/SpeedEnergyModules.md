# Conduction speed & energy modules

This adds the **speed** and **energy** readouts of the speed–precision–energy
trade-off pipeline to the Cossell & Attwell / Richardson axon model. Three new
functions plus a small change to `Model.m`.

| File | Role |
|------|------|
| `Model.m` *(modified)* | optionally exposes per-node ionic (Na⁺) currents |
| `AuxiliaryFunctions/conductionSpeed.m` | conduction velocity (m/s) |
| `AuxiliaryFunctions/energyConsumption.m` | Na⁺ charge per AP per mm (nC/mm) |
| `AuxiliaryFunctions/interiorNodeIndices.m` | shared interior-node selector |
| `measureAxon.m` | one-call pipeline returning both readouts |
| `Examples/SpeedEnergyDemo.m` | worked example on the cortical axon |

## Quick start

```matlab
par = Carcamo2017CortexAxon();
m   = measureAxon(par);

m.cv_m_per_s               % conduction velocity (m/s)
m.energy_nC_per_AP_per_mm  % Na+ charge per AP per mm (nC/mm)
m.atp_per_AP_per_mm        % ATP molecules per AP per mm
```

Or run `SpeedEnergyDemo` for a fully worked, annotated example.

## What changed in `Model.m`

The signature gained one optional input and one optional output; **existing
calls are unaffected** (default `recordCurrents = 0` reproduces the original
behaviour exactly):

```matlab
[MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR, CURRENTS] = ...
        Model(par, filename, isVerbose, recordCurrents)
```

* `recordCurrents = 0` (default): nothing recorded, no overhead.
* `recordCurrents = 1`: the time integral of each channel's current is
  accumulated per node segment and returned in `CURRENTS` (cheap).
* `recordCurrents = 2`: additionally stores the full time-resolved current
  `CURRENTS.I` (`T x N x C`, memory heavy — use only for single runs/debugging).

`CURRENTS.charge` is the `(N node-segments × C channels)` matrix of
`∫|I| dt` in nanocoulombs; `CURRENTS.channames` labels the channels
(`'Fast Na+'`, `'Persistent Na+'`, `'Slow K+'`, …).

### Why this is the right place to compute it

The Na⁺ current is `I = g_open · (V − E_Na)`. The open conductance `g_open`
depends on the gating variables, which live only inside `Model.m` and are not
saved in the voltage trace. Computing the current inside the existing
active-channel loop (reusing `tempprod` and the node voltage) is exact and
avoids re-deriving the gating from the saved voltages.

### Units

The model works in `ms, mm, mV, mS`. Hence current `= mS · mV = µA` and charge
`= µA · ms = nC`, consistent with the unit table in `simunits.m`.

## Conduction speed — `conductionSpeed.m`

Arrival time at each interior node is the first upward crossing of a threshold
(default **−20 mV**, sub-sample interpolated). Boundary nodes are discarded
(the stimulus end and the sealed far end are not at conduction steady state).

* **`regression`** (default): least-squares fit of node distance (mm) vs.
  arrival time (ms); the slope is the velocity. Robust to per-node jitter and
  reports `R²` — recommended for sweeps.
* **`twopoint`**: distance ÷ time between the first and last interior node via
  the shipped `velocities()` helper.

Because distance is in mm and time in ms, the slope mm/ms **is** m/s.
Propagation failure returns `cv = NaN` with `failed = true`.

```matlab
cv = conductionSpeed(V, IL, dt, 'Vcross', -20);          % cv.cv in m/s
cv = conductionSpeed(V, IL, dt, 'Method', 'twopoint');
```

## Energy — `energyConsumption.m`

Metabolic-cost proxy: time-integrated absolute Na⁺ current per propagating AP
per unit distance. Na⁺ channels are matched by name (substring `'Na'`, so both
fast and persistent Na⁺ are summed), integrated charge is summed over the
interior nodes, divided by the spanned axon length and by the number of APs.

```matlab
en = energyConsumption(CURRENTS, IL, 'NodeSegPerNode', par.geo.nnodeseg);
en.chargePerAPperMM   % primary readout, nC/mm
en.atpPerAPperMM      % Na+ ions / elementary charge / 3 Na+-per-ATP
```

The Na⁺ that enters during a spike is extruded by the Na⁺/K⁺-ATPase at 3 Na⁺
per ATP, so Na⁺ influx is a standard proxy for signalling energy cost
(Attwell & Laughlin, 2001). `atpPerAPperMM` makes that conversion explicit.

## Sweeping internode length (Aim 1)

`sweepInternodeLength.m` runs the deterministic L-sweep: it varies the
homogeneous internode length while holding the **number of internodes fixed**,
so total axon length floats as `N·L` and the model stays structurally identical
(same node count, same Na⁺ sites) at every L. The number of internode segments
is scaled with L so the segment *length* stays ~constant (a d_λ-style resolution
rule), preventing long internodes from being under-resolved.

Rescaling the segments rebuilds the periaxonal-space vector, which would wipe the
**paranodal seal** (the tight, high-resistance periaxonal constriction at each
internode end that drives double-cable conduction). The sweep therefore restores
the seal after every rescale (`'PreserveParanode'`, default true), holding the
paranode length constant. This matters a lot: leaving the seal off shifts the CV
peak substantially (≈80 µm without the seal vs ≈45–50 µm with it), so keep it on
unless you are deliberately studying the sealless case.

```matlab
R = sweepInternodeLength(20:10:200);          % CV + energy at each L, with a plot
[~, k] = max(R.cv_m_per_s);                    % speed optimum
fprintf('CV peaks near L = %g um\n', R.L_um(k));
```

`R` holds one row per L: `L_um`, `cv_m_per_s`, `energy_nC_per_AP_per_mm`,
`atp_per_AP_per_mm`, `r2`, `propagated`, `nIntSeg`, `totalAxonLength_mm`
(plus `R.table` if your MATLAB has `table`). Useful options:

* `'AxonFcn'` — baseline builder (e.g. `@Ford2015SBC`) to sweep a different class.
* `'ScaleSegments', false` — hold segment count fixed instead of scaling with L.
* `'MaxSegments', 200` — cap segments to speed up large-L runs (slightly coarser).
* `'Tmax_ms', 6` — lengthen the simulation if long axons fail to propagate.
* `'SaveCsv', 'sweep.csv'` — write the table to disk.

Run `SweepDemo` for a worked sweep that plots CV-vs-L and energy-vs-L with the
cortical baseline marked. **Run time grows with L** (segment count scales up), so
start with a moderate range and extend once it looks right. Any row with
`propagated = 0` means the AP did not reach all interior nodes within `tmax`;
re-run those points with a larger `'Tmax_ms'`.

## Notes & assumptions

* Like the shipped `velocities()` helper, the modules treat each column of the
  membrane-potential matrix as one node site (one segment per node,
  `nnodeseg = 1`, as in all eight bundled fiber classes). For `nnodeseg > 1`,
  pass `'NodeSegPerNode'` to `energyConsumption`; segment charges are summed
  per node.
* `conductionSpeed` and `energyConsumption` use the **same** interior-node
  selection (`interiorNodeIndices`), so speed and energy are reported over the
  same stretch of axon and remain comparable across an internode-length sweep.
* The energy readout is a Na⁺-current proxy; it approximates but does not
  directly compute ATP turnover.
