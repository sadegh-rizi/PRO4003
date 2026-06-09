# Step 2 — Heterogeneous internode lengths: implementation plan

**Goal.** Extend the step-1 clamped homogeneous pipeline to axons whose internode lengths vary
along a single fibre, at controllable levels of heterogeneity, and compare each heterogeneous
axon to the homogeneous axon of **equal mean internode length, equal total length, equal N**.
This tests H2 (does heterogeneity reach trade-off regions a homogeneous axon cannot?) and feeds
H3 (do fibre-class distributions sit near different optima?).

---

## 0. What already exists vs. what is new

| Layer | Status | File |
|---|---|---|
| Length sampler (tunable CV, clamp total, seed, base-MATLAB) | **Done** | `sampleInternodeLengths.m` |
| Per-internode length application | Exists (vector mode, case 1b) | `UpdateInternodeLength.m` |
| Deterministic CV + Na⁺ energy, geometry-agnostic | Exists | `conductionSpeed.m`, `energyConsumption.m` |
| Arrival time / threshold crossing | Exists | `arrivalTime.m` |
| Homogeneous clamped builder | Exists (scalar L) | `buildClampedAxon.m` |
| Deterministic sweep driver | Exists (sweeps L) | `sweepCVEnergy.m` |
| Noisy precision driver | Exists (sweeps L) | `sweepPrecision.m` |
| **Heterogeneous builder (vector L, clamped)** | **NEW** | `buildHeterogeneousAxon.m` |
| **Heterogeneity drivers (sweep CV_L × realizations)** | **NEW** | `sweepHetCVEnergy.m`, `sweepHetPrecision.m` |
| **Arrangement sub-experiment** | **NEW (small)** | option inside the drivers |

The measurement functions compute node positions as `pos = [0, cumsum(IL)]`, so they already
handle non-uniform internodes correctly. `Model(par,[],verbose,rec)` returns the per-internode
vector `IL`. The missing link is a builder that *installs* a length vector into `par` while
holding everything else clamped.

> **Terminology:** `CV_L` = coefficient of variation of the internode-length distribution (the
> heterogeneity knob). `CV`/`θ` = conduction velocity (a result). Don't conflate them.

---

## 1. Design decisions (the comparison protocol)

Hold the strict-isolation clamp from step 1 and add heterogeneity as the *only* new free variable.

- **Clamped (identical across every condition):** total internode span `L_total`, axon radius,
  g-ratio, node length, channel densities, all biophysical/sim parameters, and the **number of
  internodes N**.
- **Varied:** the internode-length *distribution* (its spread `CV_L`) and, as a separate knob, the
  *spatial arrangement* of a fixed length multiset.
- **Matched comparator:** for a target mean internode length `meanL`, the homogeneous baseline is
  `buildClampedAxon(par0, meanL, L_total)`, which uses `N = round(L_total/meanL)` and gives an
  actual total `N*meanL`. To match it *exactly*, generate the heterogeneous axon with
  `N` internodes and **`TotalLength = N*meanL`** (not the raw `L_total`). Then homogeneous and
  heterogeneous share N, total, and mean exactly — any difference is pure heterogeneity.

Why fix N (not just the total): the proposal's jitter argument scales σ_t ∝ √N_nodes. Fixing N
removes node count as a confound, so a change in precision is attributable to geometry, not to
having more/fewer nodes.

**Heterogeneity levels.** Sweep `CV_L ∈ {0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0}`. `CV_L = 0` reproduces
the homogeneous axon (built-in sanity check); `CV_L = 1` is the exponential/Poisson-process
("cortical-like") case; add `CV_L > 1` only if you want a callosal-like extreme. Each level is a
*distribution*, so draw **many realizations** per level and report the spread, not one axon.

---

## 2. How myelinPatterns does it — and what to copy vs. drop

Talidou & Lefebvre (`cortical.m`/`callosal.m`) and our setup differ in important ways:

- **They sample lengths that sum to a fixed myelinated length.** Cortical: throw `L_my` unit-points
  into `numMyel` bins (Poisson counts, exact integer sum). Callosal: `gamrnd(k,θ)` rounded to µm,
  *rejected until* `sum == L_my`. **Copy the idea (fixed sum); drop the mechanism.** Our
  `sampleInternodeLengths` gets the same fixed sum by continuous rescaling — no integer rounding,
  no reject-until-exact loop, no Statistics Toolbox.
- **They let N (number of sheaths) be random per axon** (`randi`), then constrain the mean. **Drop
  this** — we fix N for a clean mean-/N-matched comparison.
- **They model PMC < 100% with long bare "exposed" segments.** **Drop for the core experiment** —
  our axon is fully myelinated (internodes + ~1 µm nodes); heterogeneity lives in internode length
  only. (PMC/long-gap heterogeneity is a possible later extension; their 5%-sheath-removal
  demyelination loop is a ready template if you go there.)
- **They pad each axon with copies of the first 7 segments ×3 at both ends** to kill boundary
  artefacts, and measure only the interior "reference axon." **Copy the spirit** — but we already
  achieve this by discarding ~10% of boundary nodes (`interiorNodeIndices`) and measuring over a
  fixed physical window (`WindowFrac`). Keep using that.
- **Phenomenological model, no ionic currents → no energy axis.** This is exactly the gap we fill:
  the Richardson double-cable HH model gives Na⁺ load → ATP, so we get the third (energy) axis they
  cannot.

Net: the only myelinPatterns concept we adopt is *fixed-sum sampling*, which `sampleInternodeLengths`
already implements more cleanly.

---

## 3. Component 1 — the sampler (already built)

`sampleInternodeLengths(N, meanL, 'CV', CV_L, 'TotalLength', N*meanL, 'Seed', s)` returns a 1×N row
vector in µm that sums exactly to `N*meanL` (so mean = meanL exactly), with generating spread `CV_L`
(`CV_L=0` → all equal; `CV_L=1` → exponential/Poisson). Base-MATLAB only.

**Empirical-data alternative (preferred if you obtain it).** If you get measured internode lengths
(Call & Bergles 2021 cortical segment data; Stedehouder 2019; Chong 2012 — the set Talidou used;
Tomassy 2014), use them two ways: (a) **fit** a Gamma to estimate the real `CV_L` so your sweep
brackets the biological value; (b) **bootstrap** — resample measured lengths to N internodes and
rescale to `N*meanL`. Same builder/drivers downstream. Add a thin `sampleInternodesEmpirical.m`
later if needed; not on the critical path.

---

## 4. Component 2 — `buildHeterogeneousAxon.m` (the one new core piece)

Mirror `buildClampedAxon` but accept a **length vector** instead of a scalar.

```
function [par, info] = buildHeterogeneousAxon(par0, lengths_um, varargin)
%  lengths_um : 1xN internode lengths (µm), already summing to the clamped total.
%  Holds radius, g-ratio, node length, channels (clamp); installs per-internode lengths;
%  re-applies the paranodal seal that the node-count rebuild wipes.
```

Steps (and where they differ from the scalar builder):

1. `nIntn = numel(lengths_um); nNode = nIntn + 1;` — **N is given by the vector**, not `round(total/L)`.
2. **Segment resolution (subtlety).** `nintseg` is global, so each internode is split into the same
   number of segments and segment length = `L_i/nseg` varies. Two valid choices:
   - *Match the comparator* (recommended for fairness): `nseg = round(nseg0 * meanL / L0)` — same
     resolution the homogeneous axon at `meanL` uses.
   - *Conservative accuracy:* `nseg = round(nseg0 * max(lengths_um) / L0)` — guarantees the longest
     internode is well resolved.
   Pick one, then **verify results are stable when you double `nseg`** (§9). Reuse the
   `MinSegments`/`MaxSegments` clamps.
3. `par = UpdateNumberOfNodes(par0, nNode, 'reset', 'max', false);`
4. `par = UpdateNumberOfInternodeSegments(par, nseg);`
5. **`par = UpdateInternodeLength(par, lengths_um);`** ← the key line. Case 1b (1×#internodes vector)
   sets each internode's total length; segments within an internode stay uniform.
6. **Paranodal seal, per internode (subtlety).** `buildClampedAxon` re-seals the first/last `nPn`
   segments with `nPn = round(pnLen*nintseg/L)`. With heterogeneous lengths `nPn_i` differs per
   internode, so loop:
   ```
   for i = 1:nIntn
       nPn_i = max(1, round(pnLen_um * par.geo.nintseg / lengths_um(i)));
       nPn_i = min(nPn_i, floor(par.geo.nintseg/2));
       segs_i = [1:nPn_i, (par.geo.nintseg-nPn_i+1):par.geo.nintseg];
       par = UpdateInternodePeriaxonalSpaceWidth(par, sealW_nm, i, segs_i, 'min');
   end
   ```
   Check `UpdateInternodePeriaxonalSpaceWidth`'s signature first (read its header like we did for
   `UpdateInternodeLength`) to confirm the `(value, internodeIdx, segmentIdx, mode)` argument order
   for per-internode targeting. Fallback if per-internode sealing is awkward: use a constant `nPn`
   and accept a slightly variable physical paranode length — but per-internode is the correct,
   comparator-consistent choice.
7. Return `info` with `.nIntn, .lengths_um, .totalLength_um (=sum), .meanL_um, .CV_L_actual`.

**Acceptance test (do this first, before any sweep):** `buildHeterogeneousAxon(par0, meanL*ones(1,N))`
must produce byte-for-byte the same conduction velocity and energy as
`buildClampedAxon(par0, meanL, N*meanL)`. This is the regression test that proves the new builder is
wired correctly.

---

## 5. Component 3 — `sweepHetCVEnergy.m` (deterministic: speed + energy)

Mirror `sweepCVEnergy`, but the outer loop is over `CV_L` levels and an inner loop over realizations.
CV and energy are deterministic, so **one model run per realization** — cheap.

```
for each CV_L in levels
    for r = 1:nReal
        L = sampleInternodeLengths(N, meanL, 'CV', CV_L, 'TotalLength', N*meanL, 'Seed', baseSeed+r);
        par = buildHeterogeneousAxon(par0, L, builderOpts{:});
        [V, IL, t, C] = Model(par, [], false, 1);              % deterministic, currents on
        cv = conductionSpeed(V, IL, dt, 'NodeRange', windowFromPos(IL));
        en = energyConsumption(C, IL, 'NodeSegPerNode', par.geo.nnodeseg, 'NodeRange', sameRange);
        store cv.cv, en.chargePerAPperMM, en.atpPerAPperMM, plus realized mean/CV_L/total
    end
end
```

Reuse `sweepCVEnergy`'s fixed-physical-window logic (`WindowFrac [0.2 0.8]`) so speed and energy
come from the same interior nodes, comparable across realizations. Output one row per
(CV_L, realization); summarize as **mean ± SD across realizations** per CV_L, and keep the raw cloud
for the trade-off plot. Use the **same `baseSeed+r`** here and in the precision driver so a given
realization index refers to the same geometry in both.

Budget: e.g. 7 levels × 50 realizations = 350 fast runs. Trivial; no parallelism needed.

---

## 6. Component 4 — `sweepHetPrecision.m` (noisy: timing precision)

Mirror `sweepPrecision`. For each (CV_L, realization) build the heterogeneous axon, do one
noise-free geometry run to fix node positions and the distal node, then run `nTrials` noisy trials
(`par.noise.amp`, `par.noise.seed = trialSeedBase + k`) and take σ_t at the distal node. **Reuse
the existing noise mechanism and parfor/CSV-checkpoint machinery unchanged** — only the builder call
changes.

**Two variance sources — report both, keep them separate:**
- **Within-realization jitter** σ_t (noise only, same metric as step 1). Per CV_L, look at the
  *distribution of σ_t across realizations* and compare to the homogeneous σ_t at the same mean.
- **Across-realization spread** of the *deterministic* mean latency (geometry only; you already have
  the noise-free latency from the geometry run). This is the path-dependence Talidou & Lefebvre
  report, and it exists even with zero membrane noise.

Normalize σ_t per unit distance (`sigma_per_mm = σ_t / d`) exactly as step 1.

**Budget (this is the heavy step).** levels × realizations × trials. Keep it tractable: fewer
realizations than the deterministic sweep (e.g. 10–20) and `nTrials = 100`. 7 × 15 × 100 ≈ 10.5k
runs → run with `Parallel=true`, `SaveCsv` checkpointing, and batch by CV_L level (this maps to the
"run precision in batches" task already in your contract). **Re-tune `NoiseAmp_nA`** so jitter is
sub-millisecond at this N before the big run (§9).

---

## 7. Component 5 — arrangement sub-experiment (isolates the ordering knob)

Heterogeneity has two knobs: the length *distribution* (§5–6) and the spatial *arrangement* of those
lengths. Talidou & Lefebvre show ordering matters (velocity rises long→short, falls short→long). To
isolate it, hold one length multiset fixed and permute it:

- `random` (as sampled), `ascending`, `descending`, `clustered` (all long internodes together),
  `interleaved` (alternating long/short).

Feed each permutation through `buildHeterogeneousAxon` → CV/energy/precision. Same multiset ⇒ same
mean, total, CV_L, N; only order changes. A measurable spread across permutations confirms
path-dependence in *our* HH model and is a clean, cheap result. Implement as an `'Order'` option in
the builder or a one-line `permuteLengths(L, mode)` helper.

---

## 8. Analysis & comparison

- **Per-axis curves:** plot CV, ATP/mm, and σ_per_mm vs `CV_L` (mean ± SD across realizations), with
  the homogeneous `CV_L=0` point as reference. Tests H1 at the heterogeneity level.
- **3D trade-off cloud (the H2 test):** map each realization to the normalized point
  `m = (θ̂, 1/σ̂_t, 1/Ê)` (higher = better on all three; reuse the step-1 normalization). Plot the
  homogeneous baseline as a single point and the heterogeneous realizations as a cloud per CV_L.
  **H2 is supported if the cloud reaches regions the homogeneous point cannot** (e.g. better precision
  at equal speed, or a Pareto-nondominated zone). Quantify with the fraction of realizations that
  Pareto-dominate, or the convex-hull volume gained.
- **Direction check (Jensen):** from the step-1 CV-vs-L curve, near the optimum velocity is concave,
  so heterogeneity should on average *lower* CV at fixed mean; confirm the sign and compare with the
  arrangement spread.
- **H3 (if time):** repeat at each fibre class's empirical (mean, CV_L) — cortex ≈ low CV_L, callosum
  ≈ high CV_L (per Talidou & Lefebvre / Chong) — and see which objective each class sits nearest.

---

## 9. Verification checklist (build these as you go)

1. **Homogeneous-limit regression** (§4): `buildHeterogeneousAxon(par0, meanL*ones(1,N))` ≡
   `buildClampedAxon(par0, meanL, N*meanL)` in CV and energy. Must match to numerical tolerance.
2. **Sampler invariants** (already verified offline): `sum(L)=N*meanL`, `mean(L)=meanL`,
   realized `CV_L`→target as N grows, `CV_L=0` ⇒ all equal. Add a `assert` in the driver.
3. **Resolution convergence:** rerun one CV_L level with `nseg` doubled; CV/energy should be stable
   (<~1%). If not, raise `nseg` (use the max-length rule).
4. **Noise recalibration:** with the step-2 N, sweep `NoiseAmp_nA` on the homogeneous axon until σ_t
   is sub-ms and matches your step-1 jitter regime, *before* the big precision run.
5. **Trend sanity vs. literature:** more heterogeneity ⇒ larger dispersion in delay/CV and (at low
   coverage / long internodes) more failures — consistent with Talidou & Lefebvre.
6. **Reproducibility:** same `Seed` ⇒ identical geometry; record seeds in the output CSV.

---

## 10. Sequencing & task mapping (fits your group contract)

1. `buildHeterogeneousAxon.m` + regression test (§4, §9.1) — **coding pair (Sadegh/Pablo)**. Blocks everything.
2. `sweepHetCVEnergy.m` + run (§5) — deterministic, fast — **Sadegh** (CV/energy module owner).
3. Noise recalibration (§9.4) + `sweepHetPrecision.m` (§6) — **Tess/Pablo** (precision + noise owners).
4. Arrangement helper + sub-experiment (§7) — small, **Pablo/Sadegh**.
5. Empirical-data fit/bootstrap (§3) — **Tess/Ioannis** (data search already in the contract).
6. Analysis, 3D trade-off, H2/H3 figures (§8) — **Pablo (+Sadegh)**, mirrors the "combined objective" task.

Critical path: 1 → 2 → 3. Items 4–5 run in parallel once 1 is done.

## 11. Pitfalls

- Don't compare heterogeneous to a homogeneous axon at the *raw* total — use `N*meanL` so N, total,
  and mean all match (§1).
- A global `nintseg` makes long internodes coarser; confirm convergence (§9.3).
- Forgetting to re-seal paranodes after the node-count rebuild changes conduction — keep
  `PreserveParanode` behaviour (§4.6).
- `gamrnd` (myelinPatterns) needs the Statistics Toolbox; our base-MATLAB sampler does not — don't
  reintroduce the dependency.
- Precision is the compute sink: batch by CV_L, checkpoint with `SaveCsv`, re-tune noise first.

---

### References
- Talidou & Lefebvre 2025, *eNeuro*, doi:10.1523/ENEURO.0402-24.2025; code: github.com/atalidou/myelinPatterns (`cortical.m` = Poisson, `callosal.m` = Gamma).
- Internode-length data to fit/bootstrap: Chong et al. 2012 *PNAS*; Call & Bergles 2021 *Nat Commun*; Stedehouder et al. 2019 *eLife*; Tomassy et al. 2014 *Science*.
