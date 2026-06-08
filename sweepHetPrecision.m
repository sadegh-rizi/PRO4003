function R = sweepHetPrecision(meanL_um, CVlevels, varargin)
%SWEEPHETPRECISION  Timing precision vs heterogeneity level (fixed mean, clamped total, noisy trials).
%
%   R = SWEEPHETPRECISION(meanL_um, CVlevels, ...)
%
%   Step-2 analogue of sweepPrecision. The mean internode length is held fixed
%   and the spread (CV_L) is swept. For each CV_L level it draws nReal
%   heterogeneous axons (sum clamped to N*meanL, N fixed) and, for each, runs
%   one noise-free geometry pass (to fix node positions and the distal node)
%   then nTrials noisy trials. It reports TWO distinct variabilities:
%       * within-realization jitter  sigma_t  = std of distal arrival time over
%         the noisy trials (the step-1 precision metric), and
%       * across-realization geometry latency detLatency_ms = the noise-free
%         distal arrival time (path-dependence; varies between draws even with
%         no membrane noise).
%   N is fixed, so the proposal's sigma_t ~ sqrt(N_nodes) scaling is controlled
%   and differences reflect geometry, not node count. CV_L = 0 is the homogeneous
%   baseline.
%
%   Name/value options:
%       'nReal'          - heterogeneous realizations per CV_L level (default 15).
%       'nTrials'        - noisy trials per realization (default 100).
%       'NoiseAmp_nA'    - noise intensity (default 0.05). RETUNE for sub-ms jitter
%                          at your N before a production run (see plan §9.4).
%       'TotalLength_um' - clamped total internode span (default = baseline N*L0).
%       'AxonFcn'        - baseline axon builder (default @Carcamo2017CortexAxon).
%       'Order'          - arrangement of each draw (see permuteLengths; default 'random').
%       'SamplerMin'/'SamplerMax' - length bounds for the sampler (default 0/Inf).
%       'BaseSeed'       - geometry/sampler seed = BaseSeed + 1000*levelIdx + real (default 7000).
%                          (Match sweepHetCVEnergy's BaseSeed/CVlevels for identical geometries.)
%       'TrialSeedBase'  - noise seed = TrialSeedBase + 100000*levelIdx + 1000*real + trial (default 1e6).
%       'Vcross'         - arrival threshold, mV (default -20).
%       'DistalFrac'     - distal node position as a fraction of total (default 0.8).
%       'ExcludeStimNode'- keep noise off the stimulated node (default false).
%       'Tmax_ms'/'Dt_us'- sim duration / time step overrides.
%       'Parallel'       - run trials over parfor (default true; serial if no pool).
%       'Verbose'        - print a line per realization (default true).
%       'SaveCsv'        - path to the flat results table (default []).
%       (builder options Nseg/ResolveBy/MinSegments/MaxSegments/... forwarded.)
%
%   Output struct R:
%       .CV_L, .nReal, .nTrials, .noiseAmp_nA, .N, .meanL_um, .totalLength_um,
%       .sigma_t (nLevels x nReal, ms), .sigma_per_mm (ms/mm),
%       .detLatency_ms, .meanLatency_ms, .distance_mm, .nValid,
%       per-level .sigma_t_mean/.sigma_t_std, .sigma_per_mm_mean,
%       .detLatency_std (across-realization geometry spread), plus .table.
%
%   NOTE: heavy step (levels x nReal x nTrials model runs). Batch by level and
%   use 'SaveCsv'. See also: sweepPrecision, buildHeterogeneousAxon, arrivalTime.

axonFcn   = getOption(varargin, 'AxonFcn', @Carcamo2017CortexAxon);
nReal     = getOption(varargin, 'nReal', 15);
nTrials   = getOption(varargin, 'nTrials', 100);
noiseAmp  = getOption(varargin, 'NoiseAmp_nA', 0.05);
totalLen  = getOption(varargin, 'TotalLength_um', []);
order     = getOption(varargin, 'Order', 'random');
sMin      = getOption(varargin, 'SamplerMin', 0);
sMax      = getOption(varargin, 'SamplerMax', inf);
baseSeed  = getOption(varargin, 'BaseSeed', 7000);
trialBase = getOption(varargin, 'TrialSeedBase', 1e6);
vcross    = getOption(varargin, 'Vcross', -20);
distFrac  = getOption(varargin, 'DistalFrac', 0.8);
exclStim  = getOption(varargin, 'ExcludeStimNode', false);
tmaxMs    = getOption(varargin, 'Tmax_ms', []);
dtUs      = getOption(varargin, 'Dt_us', []);
parWanted = getOption(varargin, 'Parallel', true);
verbose   = getOption(varargin, 'Verbose', true);
saveCsv   = getOption(varargin, 'SaveCsv', []);
builderOpts = builderOptions(varargin);

Mworkers = Inf; if ~parWanted, Mworkers = 0; end

par0 = axonFcn();
if isempty(totalLen)
    totalLen = par0.geo.nintn * par0.intn.geo.length.value.ref;
end
N           = max(1, round(totalLen / meanL_um));
targetTotal = N * meanL_um;

CVlevels = CVlevels(:).';
nLev = numel(CVlevels);

sig   = nan(nLev, nReal);  signorm = nan(nLev, nReal);
detL  = nan(nLev, nReal);  meanL   = nan(nLev, nReal);
dist  = nan(nLev, nReal);  nValid  = zeros(nLev, nReal);

if verbose
    fprintf('Het precision sweep: meanL=%.1f um, N=%d, %d levels x %d real x %d trials, noise=%g nA\n', ...
            meanL_um, N, nLev, nReal, nTrials, noiseAmp);
    fprintf('%6s %5s %9s %9s %12s %14s %6s\n', 'CV_L', 'real', 'CV_L(act)', 'd[mm]', 'sigma_t[ms]', 'sigma[ms/mm]', 'nVal');
end

for ci = 1:nLev
    CVL = CVlevels(ci);
    for r = 1:nReal
        try
            seed = baseSeed + 1000 * ci + r;
            L = sampleInternodeLengths(N, meanL_um, 'CV', CVL, 'TotalLength', targetTotal, ...
                                       'Min', sMin, 'Max', sMax, 'Seed', seed);
            if ~strcmpi(order, 'random'), L = permuteLengths(L, order); end

            par = buildHeterogeneousAxon(par0, L, builderOpts{:});
            if ~isempty(tmaxMs), par.sim.tmax.value = tmaxMs; par.sim.tmax.units = {1,'ms',1}; end
            if ~isempty(dtUs),   par.sim.dt.value   = dtUs;   par.sim.dt.units   = {1,'us',1};  end

            % Geometry pass (noise OFF): node positions, distal node, det. latency.
            [Vg, IL, tg] = Model(par, [], false, 0);
            dt  = tg(2) - tg(1);
            pos = [0, cumsum(IL(:).')];
            [~, distalNode] = min(abs(pos - distFrac * pos(end)));
            distalNode = min(max(distalNode, 2), numel(pos) - 1);
            d = pos(distalNode);
            detLat = arrivalTime(Vg(:, distalNode), dt, vcross);

            % Noisy trials (independent, reproducible).
            arrivals = nan(1, nTrials);
            seedT0 = trialBase + 100000 * ci + 1000 * r;
            parfor (k = 1:nTrials, Mworkers)
                park = par;
                park.noise.amp.value      = noiseAmp;
                park.noise.amp.units      = {1, 'nA', 1};
                park.noise.excludeStimNode = exclStim;
                park.noise.seed            = seedT0 + k;
                Vk = Model(park, [], false, 0);
                arrivals(k) = arrivalTime(Vk(:, distalNode), dt, vcross);
            end

            good = ~isnan(arrivals);
            dist(ci, r)   = d;       detL(ci, r) = detLat;
            nValid(ci, r) = sum(good);
            meanL(ci, r)  = mean(arrivals(good));
            sig(ci, r)    = std(arrivals(good));
            signorm(ci, r)= sig(ci, r) / d;

            if verbose
                fprintf('%6.2f %5d %9.3f %9.3f %12.4g %14.4g %6d\n', ...
                        CVL, r, std(L)/mean(L), d, sig(ci, r), signorm(ci, r), nValid(ci, r));
            end
        catch err
            if verbose, fprintf('  CV_L=%.2f real=%d FAILED: %s\n', CVL, r, err.message); end
        end
    end
    if ~isempty(saveCsv), writeFlat(saveCsv, CVlevels, sig, signorm, detL, meanL, dist, nValid); end
end

R = struct();
R.CV_L = CVlevels(:);  R.nReal = nReal;  R.nTrials = nTrials;  R.noiseAmp_nA = noiseAmp;
R.N = N;  R.meanL_um = meanL_um;  R.totalLength_um = targetTotal;
R.sigma_t = sig;  R.sigma_per_mm = signorm;  R.detLatency_ms = detL;
R.meanLatency_ms = meanL;  R.distance_mm = dist;  R.nValid = nValid;
R.sigma_t_mean     = nanmeanRow(sig);      R.sigma_t_std      = nanstdRow(sig);
R.sigma_per_mm_mean = nanmeanRow(signorm); R.detLatency_std   = nanstdRow(detL);

[CIg, Rg] = ndgrid(1:nLev, 1:nReal);
try
    R.table = table(reshape(CVlevels(CIg(:)), [], 1), Rg(:), sig(:), signorm(:), detL(:), meanL(:), dist(:), nValid(:), ...
        'VariableNames', {'CV_L','realization','sigma_t_ms','sigma_per_mm','detLatency_ms', ...
                          'meanLatency_ms','distance_mm','nValid'});
catch
    R.table = [];
end
if ~isempty(saveCsv) && ~isempty(R.table)
    try, writetable(R.table, saveCsv); catch e, warning('sweepHetPrecision:csv','%s', e.message); end
end
end


% ===================== local helpers =====================
function writeFlat(fname, CVlevels, sig, signorm, detL, meanL, dist, nValid)
fid = fopen(fname, 'w'); if fid < 0, return, end
fprintf(fid, 'CV_L,realization,sigma_t_ms,sigma_per_mm,detLatency_ms,meanLatency_ms,distance_mm,nValid\n');
for ci = 1:numel(CVlevels)
    for r = 1:size(sig,2)
        fprintf(fid, '%g,%d,%g,%g,%g,%g,%g,%d\n', CVlevels(ci), r, sig(ci,r), signorm(ci,r), ...
                detL(ci,r), meanL(ci,r), dist(ci,r), nValid(ci,r));
    end
end
fclose(fid);
end

function m = nanmeanRow(X)
m = nan(size(X,1),1);
for i = 1:size(X,1), v = X(i,~isnan(X(i,:))); if ~isempty(v), m(i) = mean(v); end, end
end

function s = nanstdRow(X)
s = nan(size(X,1),1);
for i = 1:size(X,1), v = X(i,~isnan(X(i,:))); if numel(v) > 1, s(i) = std(v); end, end
end

function opts = builderOptions(args)
names = {'Nseg','ResolveBy','MinSegments','MaxSegments','PreserveParanode', ...
         'ParanodeLength_um','ParanodeSealWidth_nm'};
opts = {};
for k = 1:2:numel(args) - 1
    if any(strcmpi(args{k}, names)), opts(end+1:end+2) = {args{k}, args{k+1}}; end %#ok<AGROW>
end
end

function val = getOption(args, name, default)
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name), val = args{k + 1}; end
end
end
