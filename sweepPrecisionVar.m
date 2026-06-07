function R = sweepPrecision2(Lvalues_um, varargin)
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
%       var_t         = var of distal arrival time (ms^2)
%       var_per_mm    = var_t / d             (ms^2/mm)
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
%       'Parallel'        - run the per-length trials over a parfor (default true).
%       'Progress'        - show an in-place per-length progress bar with elapsed
%                           time and ETA (default true). Set false for quiet logs.
%       'nBoot'           - bootstrap resamples for the 95% CI on sigma_t (default 2000).

axonFcn   = getOption(varargin, 'AxonFcn', @Carcamo2017CortexAxon);
totalLen  = getOption(varargin, 'TotalLength_um', []);
noiseAmp  = getOption(varargin, 'NoiseAmp_nA', 0.02);
nTrials   = getOption(varargin, 'nTrials', 2);
vcross    = getOption(varargin, 'Vcross', -20);
distFrac  = getOption(varargin, 'DistalFrac', 0.8);
exclStim  = getOption(varargin, 'ExcludeStimNode', false);
baseSeed  = getOption(varargin, 'BaseSeed', 1000);
verbose   = getOption(varargin, 'Verbose', true);
saveCsv   = getOption(varargin, 'SaveCsv', []);
tmaxMs    = getOption(varargin, 'Tmax_ms', []);
dtUs      = getOption(varargin, 'Dt_us', []);
parWanted = getOption(varargin, 'Parallel', true);
showProg  = getOption(varargin, 'Progress', true);
nBoot     = getOption(varargin, 'nBoot', 2000);     

Mworkers = Inf;
if ~parWanted, Mworkers = 0; end

par0 = axonFcn();
if isempty(totalLen)
    totalLen = par0.geo.nintn * par0.intn.geo.length.value.ref;
end

Lvalues_um = Lvalues_um(:).';
nL = numel(Lvalues_um);

L_act = nan(nL,1); nIntn = nan(nL,1); nNode = nan(nL,1); dNode = nan(nL,1);
dist  = nan(nL,1); meanLat = nan(nL,1); nValid = zeros(nL,1); failed_rate = nan(nL,1);   

% Arrays for Standard Deviation and Variance
sig = nan(nL,1); signorm = nan(nL,1); 
var_t = nan(nL,1); varnorm = nan(nL,1);

sigLo = nan(nL,1); sigHi = nan(nL,1);   
snLo  = nan(nL,1); snHi  = nan(nL,1);   

if verbose
    fprintf('Precision sweep: total = %.0f um, %d lengths x %d trials, noise = %g nA\n', ...
            totalLen, nL, nTrials, noiseAmp);
    if ~showProg
        fprintf('%7s %6s %9s %12s %14s %14s\n','L[um]','nIntn','d[mm]','sigma_t[ms]','sigma[ms/mm]','var[ms2/mm]');
    end
end
t0 = tic;   

for q = 1:nL
    L = Lvalues_um(q);
    try
        % [par, info] = buildClampedAxon(par0, L, totalLen, builderOpts{:});
        [par, info] = buildClampedAxon(par0, L, totalLen);
        if ~isempty(tmaxMs), par.sim.tmax.value = tmaxMs; par.sim.tmax.units = {1,'ms',1}; end
        if ~isempty(dtUs),   par.sim.dt.value   = dtUs;   par.sim.dt.units   = {1,'us',1};  end
        par.noise.amp.value      = noiseAmp;
        par.noise.amp.units      = {1, 'nA', 1};
        par.noise.excludeStimNode = exclStim;

        % Spatial Geometry Interpolation Setup
        parGeom = par;  parGeom.noise.amp.value = 0;        
        [~, IL, tg] = Model(parGeom, [], false, 0);
        dt  = tg(2) - tg(1);
        pos = [0, cumsum(IL(:).')];                         
        
        dTarget = distFrac * pos(end);                      
        iR = find(pos >= dTarget, 1, 'first');              
        iR = min(max(iR, 2), numel(pos));
        iL = iR - 1;                                        
        
        if pos(iR) == pos(iL)
            fracPos = 0; 
        else
            fracPos = (dTarget - pos(iL)) / (pos(iR) - pos(iL)); 
        end
        
        d = dTarget;                                        
        distalNode = iL + fracPos;                          

        arrivals = nan(1, nTrials);
        if showProg
            nChunks = max(1, min(10, floor(nTrials / 5)));
            fprintf('[%d/%d] L=%5.1f um  %3.0f%%', q, nL, L, 0);
        else
            nChunks = 1;
        end
        edges = round(linspace(0, nTrials, nChunks + 1));
        
        for c = 1:nChunks
            idx    = (edges(c) + 1):edges(c + 1);
            aChunk = nan(1, numel(idx));
            
            parfor (j = 1:numel(idx), Mworkers)
                kk = idx(j);
                park = par;
                park.noise.seed = baseSeed + kk;            
                Vk = Model(park, [], false, 0);             
                
                tL = arrivalTime(Vk(:, iL), dt, vcross);
                tR = arrivalTime(Vk(:, iR), dt, vcross);
                
                if ~isnan(tL) && ~isnan(tR)
                    aChunk(j) = tL + fracPos * (tR - tL);
                else
                    aChunk(j) = NaN; 
                end
            end
            
            arrivals(idx) = aChunk;
            if showProg, fprintf('\b\b\b\b%3.0f%%', 100 * edges(c + 1) / nTrials); end
        end
        
        good = ~isnan(arrivals);
        L_act(q)  = info.L_um;  nIntn(q) = info.nIntn;  nNode(q) = info.nNode;
        dNode(q)  = distalNode; dist(q)  = d;            nValid(q) = sum(good);
        meanLat(q)= mean(arrivals(good));
        
        % Calculate BOTH Standard Deviation and Variance
        sig(q)     = std(arrivals(good));
        var_t(q)   = var(arrivals(good)); % or sig(q)^2
        
        % Normalize both by distance
        signorm(q) = sig(q) / d;
        varnorm(q) = var_t(q) / d;
        
        failed_rate(q) = 1 - nValid(q) / nTrials;

        a = arrivals(good);
        if numel(a) >= 3
            rng(baseSeed + 100000 + q);             
            bs = zeros(1, nBoot);
            for b = 1:nBoot
                bs(b) = std(a(randi(numel(a), 1, numel(a))));
            end
            bs = sort(bs);
            sigLo(q) = bs(max(1,     round(0.025 * nBoot)));
            sigHi(q) = bs(min(nBoot, round(0.975 * nBoot)));
            snLo(q)  = sigLo(q) / d;
            snHi(q)  = sigHi(q) / d;
        end
        
        if showProg   
            el = toc(t0);  eta = el / q * (nL - q);
            fprintf('\b\b\b\b  sigma_t=%.3g ms  var/mm=%.3g  [elapsed %.1f min, ETA %.1f min]\n', ...
                    sig(q), varnorm(q), el/60, eta/60);
        end
        
    catch err
        if showProg, fprintf('\n'); end
        if verbose || showProg, fprintf('  L = %g um FAILED: %s\n', L, err.message); end
    end
    
    if verbose && ~showProg
        fprintf('%7.1f %6d %9.3f %12.4g %14.4g %14.4g\n', L, nIntn(q), dist(q), sig(q), signorm(q), varnorm(q));
    end
    
    if ~isempty(saveCsv)
        writePrecisionCsv(saveCsv, Lvalues_um, L_act, nIntn, nNode, dNode, dist, ...
                          meanLat, sig, signorm, var_t, varnorm, nValid, failed_rate, snLo, snHi);
    end
end

R = struct();
R.L_um = Lvalues_um(:); R.L_actual_um = L_act; R.nIntn = nIntn; R.nNode = nNode;
R.distalNode = dNode; R.distance_mm = dist; R.meanLatency_ms = meanLat;
R.sigma_t_ms = sig; R.sigma_per_mm = signorm; 
R.var_t_ms2 = var_t; R.var_per_mm = varnorm; % Save Variance to Output Struct
R.nValid = nValid; R.failed_rate = failed_rate;
R.sigma_t_lo = sigLo; R.sigma_t_hi = sigHi;             
R.sigma_per_mm_lo = snLo; R.sigma_per_mm_hi = snHi;     
R.clampedTotal_um = totalLen; R.nTrials = nTrials; R.noiseAmp_nA = noiseAmp;

try
    R.table = table(R.L_um, R.L_actual_um, R.nIntn, R.distance_mm, ...
                    R.meanLatency_ms, R.sigma_t_ms, R.sigma_per_mm, ...
                    R.var_t_ms2, R.var_per_mm, ... % Save Variance to Output Table
                    R.sigma_per_mm_lo, R.sigma_per_mm_hi, R.nValid, R.failed_rate, ...
        'VariableNames', {'L_um','L_actual_um','nIntn','distance_mm', ...
                          'meanLatency_ms','sigma_t_ms','sigma_per_mm', ...
                          'var_t_ms2','var_per_mm',...
                          'sigma_per_mm_lo','sigma_per_mm_hi','nValid','failed_rate'});
catch
    R.table = [];
end
end

% ===================== local helpers =====================
function writePrecisionCsv(fname, L, La, nI, nN, dN, d, ml, s, sn, v, vn, nv, failed_rate, snLo, snHi)
fid = fopen(fname, 'w');
if fid < 0, return, end
% Updated CSV Header
fprintf(fid, 'L_um,L_actual_um,nIntn,nNode,distalNode,distance_mm,meanLatency_ms,sigma_t_ms,sigma_per_mm,var_t_ms2,var_per_mm,sigma_per_mm_lo,sigma_per_mm_hi,nValid,failed_rate\n');
for i = 1:numel(L)
    % Updated CSV Data Writer
    fprintf(fid, '%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n', ...
            L(i), La(i), nI(i), nN(i), dN(i), d(i), ml(i), s(i), sn(i), v(i), vn(i), snLo(i), snHi(i), nv(i), failed_rate(i));
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