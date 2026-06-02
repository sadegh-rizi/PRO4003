function R = sweepPrecision(Lvalues_um, varargin)
%SWEEPPRECISION  Timing precision across internode lengths (clamped total, noisy trials).
%
%   R = SWEEPPRECISION(Lvalues_um, ...)
%
%   For each internode length L it builds the SAME clamped-total-length axon as
%   sweepCVEnergy (see buildClampedAxon), then runs nTrials noisy simulations
%   (independent Gaussian node-current noise via par.noise) and measures the
%   distal node's arrival time on each trial using the shared arrivalTime()
%   threshold crossing. Timing precision is the across-trial standard deviation
%   of that arrival time, normalised per unit propagation distance:
%       sigma_t       = std of distal arrival time (ms)
%       sigma_per_mm  = sigma_t / d           (ms/mm),  d = distance to distal node.
%
%   Name/value options:
%       'AxonFcn'         - baseline axon builder (default @Carcamo2017CortexAxon).
%       'TotalLength_um'  - clamped total internode span. Default = baseline N*L0.
%       'NoiseAmp_nA'     - noise intensity (default 0.05). MUST be tuned to give
%                           sub-ms jitter; note it is a dt-invariant intensity.
%       'nTrials'         - trials per length (default 100).
%       'Vcross'          - arrival threshold, mV (default -20).
%       'DistalFrac'      - distal node position as a fraction of total (default 0.8).
%       'ExcludeStimNode' - keep noise off the stimulated node (default false).
%       'BaseSeed'        - trial seeds are BaseSeed + trial index (default 1000).
%       'Verbose'         - print a line per L (default true).
%       'SaveCsv'         - path to a CSV; rewritten after each L (crash-safe batching).
%       'Tmax_ms'         - override sim duration (ms); raise if the distal node never
%                           fires (nValid = 0). Default keeps the axon's own tmax.
%       'Dt_us'           - override the time step (us); coarser = faster (re-tune noise).
%       (builder options ScaleSegments/MinSegments/MaxSegments/PreserveParanode
%        are forwarded to buildClampedAxon.)
%
%   Output struct R (one row per L): .L_um, .L_actual_um, .nIntn, .nNode,
%       .distalNode, .distance_mm, .meanLatency_ms, .sigma_t_ms,
%       .sigma_per_mm (ms/mm), .nValid, plus .table.
%
%   NOTE: this is the heavy step (nTrials x nL model runs). Run in batches by
%   calling with subsets of Lvalues_um; with 'SaveCsv' set, partial results are
%   written after every length.
%
%   See also: buildClampedAxon, arrivalTime, sweepCVEnergy, Model.

axonFcn   = getOption(varargin, 'AxonFcn', @Carcamo2017CortexAxon);
totalLen  = getOption(varargin, 'TotalLength_um', []);
noiseAmp  = getOption(varargin, 'NoiseAmp_nA', 0.05);
nTrials   = getOption(varargin, 'nTrials', 100);
vcross    = getOption(varargin, 'Vcross', -20);
distFrac  = getOption(varargin, 'DistalFrac', 0.8);
exclStim  = getOption(varargin, 'ExcludeStimNode', false);
baseSeed  = getOption(varargin, 'BaseSeed', 1000);
verbose   = getOption(varargin, 'Verbose', true);
saveCsv   = getOption(varargin, 'SaveCsv', []);
tmaxMs    = getOption(varargin, 'Tmax_ms', []);
dtUs      = getOption(varargin, 'Dt_us', []);
builderOpts = builderOptions(varargin);

par0 = axonFcn();
if isempty(totalLen)
    totalLen = par0.geo.nintn * par0.intn.geo.length.value.ref;
end

Lvalues_um = Lvalues_um(:).';
nL = numel(Lvalues_um);

L_act = nan(nL,1); nIntn = nan(nL,1); nNode = nan(nL,1); dNode = nan(nL,1);
dist  = nan(nL,1); meanLat = nan(nL,1); sig = nan(nL,1); signorm = nan(nL,1); nValid = zeros(nL,1);

if verbose
    fprintf('Precision sweep: total = %.0f um, %d lengths x %d trials, noise = %g nA\n', ...
            totalLen, nL, nTrials, noiseAmp);
    fprintf('%7s %6s %9s %12s %14s\n','L[um]','nIntn','d[mm]','sigma_t[ms]','sigma[ms/mm]');
end

for q = 1:nL
    L = Lvalues_um(q);
    try
        [par, info] = buildClampedAxon(par0, L, totalLen, builderOpts{:});
        if ~isempty(tmaxMs), par.sim.tmax.value = tmaxMs; par.sim.tmax.units = {1,'ms',1}; end
        if ~isempty(dtUs),   par.sim.dt.value   = dtUs;   par.sim.dt.units   = {1,'us',1};  end
        par.noise.amp.value      = noiseAmp;
        par.noise.amp.units      = {1, 'nA', 1};
        par.noise.excludeStimNode = exclStim;

        arrivals    = nan(1, nTrials);
        distalNode  = NaN; d = NaN; dt = NaN;
        for k = 1:nTrials
            par.noise.seed = baseSeed + k;                 % independent, reproducible
            [V, IL, t] = Model(par, [], false, 0);          % noisy, currents off
            if k == 1
                dt  = t(2) - t(1);
                pos = [0, cumsum(IL(:).')];                 % mm, one per node
                Ltot_actual = pos(end);
                [~, distalNode] = min(abs(pos - distFrac * Ltot_actual));
                distalNode = min(max(distalNode, 2), numel(pos) - 1);  % keep interior
                d = pos(distalNode);                        % mm from stimulated node 1
            end
            arrivals(k) = arrivalTime(V(:, distalNode), dt, vcross);
        end

        good = ~isnan(arrivals);
        L_act(q)  = info.L_um;  nIntn(q) = info.nIntn;  nNode(q) = info.nNode;
        dNode(q)  = distalNode; dist(q)  = d;            nValid(q) = sum(good);
        meanLat(q)= mean(arrivals(good));
        sig(q)    = std(arrivals(good));
        signorm(q)= sig(q) / d;
        fail_rate = 1 - nValid(q) / nTrials;
        
    catch err
        if verbose, fprintf('  L = %g um FAILED: %s\n', L, err.message); end
    end
    if verbose
        fprintf('%7.1f %6d %9.3f %12.4g %14.4g\n', L, nIntn(q), dist(q), sig(q), signorm(q));
    end

    % Crash-safe checkpoint after each length.
    if ~isempty(saveCsv)
        writePrecisionCsv(saveCsv, Lvalues_um, L_act, nIntn, nNode, dNode, dist, ...
                          meanLat, sig, signorm, nValid);
    end
end

R = struct();
R.L_um = Lvalues_um(:); R.L_actual_um = L_act; R.nIntn = nIntn; R.nNode = nNode;
R.distalNode = dNode; R.distance_mm = dist; R.meanLatency_ms = meanLat;
R.sigma_t_ms = sig; R.sigma_per_mm = signorm; R.nValid = nValid;
R.clampedTotal_um = totalLen; R.nTrials = nTrials; R.noiseAmp_nA = noiseAmp;
try
    R.table = table(R.L_um, R.L_actual_um, R.nIntn, R.distance_mm, ...
                    R.meanLatency_ms, R.sigma_t_ms, R.sigma_per_mm, R.nValid, ...
        'VariableNames', {'L_um','L_actual_um','nIntn','distance_mm', ...
                          'meanLatency_ms','sigma_t_ms','sigma_per_mm','nValid'});
catch
    R.table = [];
end
end


% ===================== local helpers =====================

function writePrecisionCsv(fname, L, La, nI, nN, dN, d, ml, s, sn, nv)
fid = fopen(fname, 'w');
if fid < 0, return, end
fprintf(fid, 'L_um,L_actual_um,nIntn,nNode,distalNode,distance_mm,meanLatency_ms,sigma_t_ms,sigma_per_mm,nValid\n');
for i = 1:numel(L)
    fprintf(fid, '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n', ...
            L(i), La(i), nI(i), nN(i), dN(i), d(i), ml(i), s(i), sn(i), nv(i));
end
fclose(fid);
end


function opts = builderOptions(args)
names = {'ScaleSegments','MinSegments','MaxSegments','PreserveParanode', ...
         'ParanodeLength_um','ParanodeSealWidth_nm'};
opts = {};
for k = 1:2:numel(args) - 1
    if any(strcmpi(args{k}, names))
        opts(end+1:end+2) = {args{k}, args{k+1}}; %#ok<AGROW>
    end
end
end


function val = getOption(args, name, default)
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name)
        val = args{k + 1};
    end
end
end
