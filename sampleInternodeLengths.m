function lengths = sampleInternodeLengths(nInternodes, meanLength, varargin)
%SAMPLEINTERNODELENGTHS  Heterogeneous internode lengths with controllable spread.
%   LENGTHS = SAMPLEINTERNODELENGTHS(NINTERNODES, MEANLENGTH) returns a
%   1 x NINTERNODES row vector of positive internode lengths drawn from a
%   Gamma distribution with the requested mean. By default CV = 1, which is
%   the EXPONENTIAL (Poisson-process) case -- identical in distribution to
%   the original samplePoisson.m. Set 'CV' to tune the level of heterogeneity.
%
%   Why Gamma instead of a pure exponential?
%     * A pure exponential has a FIXED coefficient of variation, CV = 1, so it
%       cannot represent "different levels of heterogeneity" (the proposal's
%       step 2). Gamma fixes the mean and lets you sweep CV via the shape k:
%           k = 1/CV^2 ,  theta = mean/k ,  CV = 1/sqrt(k).
%       CV -> 0 approaches a homogeneous axon; CV = 1 recovers the exponential
%       (maximal, Poisson-process heterogeneity). Biological internode
%       distributions are typically under-dispersed (CV < 1).
%
%   Name/value options:
%     'CV'          Coefficient of variation (std/mean) of the lengths.
%                   Default 1 (exponential / Poisson process). CV = 0 returns
%                   an exactly homogeneous vector (mean in every entry).
%     'TotalLength' If set, the sampled vector is rescaled so SUM(lengths)
%                   equals this value EXACTLY. Use this to honour the
%                   strict-isolation total-length clamp: it guarantees the
%                   heterogeneous axon has the same total span AND the same
%                   mean internode length (= TotalLength/nInternodes) as its
%                   homogeneous counterpart, removing the truncation bias that
%                   afflicts rejection sampling. Default [] (no rescale).
%     'Min'         Lower bound (um). Out-of-range draws are resampled.
%                   Default 0 (no lower bound). NOTE: bounds bias both the mean
%                   and the CV -- prefer 'CV' + 'TotalLength' over tight bounds.
%     'Max'         Upper bound (um). Default Inf.
%     'Seed'        Integer seed for reproducibility. The global RNG state is
%                   saved and restored, so calling with a seed has no lasting
%                   side effect. Default [] (use current global RNG).
%     'MaxAttempts' Safety cap on resampling tries per internode when bounds
%                   are active, to avoid an infinite loop if [Min,Max] is in a
%                   near-zero-probability region. Default 10000.
%
%   Output:
%     lengths       1 x NINTERNODES row vector (um).
%
%   Examples:
%     % Exponential/Poisson case, exactly as the old samplePoisson(50,60):
%     L = sampleInternodeLengths(50, 60);
%
%     % Mild heterogeneity, total length clamped to 3000 um (mean = 60 um):
%     L = sampleInternodeLengths(50, 60, 'CV', 0.3, 'TotalLength', 3000);
%
%     % Reproducible draw:
%     L = sampleInternodeLengths(50, 60, 'CV', 0.5, 'Seed', 42);
%
%   Feed the result straight into the model with a per-internode vector:
%     par = UpdateInternodeLength(par, L);   % case 1b: 1 x #internodes
%
%   See also UPDATEINTERNODELENGTH, BUILDCLAMPEDAXON.

    narginchk(2, inf);

    % ----- options -----
    CV          = getOption(varargin, 'CV',          1);
    totalLength = getOption(varargin, 'TotalLength', []);
    minLength   = getOption(varargin, 'Min',         0);
    maxLength   = getOption(varargin, 'Max',         inf);
    seed        = getOption(varargin, 'Seed',        []);
    maxAttempts = getOption(varargin, 'MaxAttempts', 10000);

    % ----- validation -----
    if ~isscalar(nInternodes) || nInternodes <= 0 || floor(nInternodes) ~= nInternodes
        error('sampleInternodeLengths:nInternodes', 'nInternodes must be a positive integer.');
    end
    if ~isscalar(meanLength) || meanLength <= 0
        error('sampleInternodeLengths:meanLength', 'meanLength must be a positive scalar.');
    end
    if ~isscalar(CV) || CV < 0
        error('sampleInternodeLengths:CV', 'CV must be a non-negative scalar.');
    end
    if minLength < 0
        error('sampleInternodeLengths:Min', 'Min must be >= 0.');
    end
    if maxLength < minLength
        error('sampleInternodeLengths:Max', 'Max must be >= Min.');
    end
    if ~isempty(totalLength) && (~isscalar(totalLength) || totalLength <= 0)
        error('sampleInternodeLengths:TotalLength', 'TotalLength must be a positive scalar.');
    end

    % ----- reproducibility (restore global RNG on exit) -----
    if ~isempty(seed)
        prevState = rng;                                 %#ok<NASGU>
        restore   = onCleanup(@() rng(prevState));        %#ok<NASGU>
        rng(seed);
    end

    % ----- homogeneous limit -----
    if CV == 0
        lengths = meanLength * ones(1, nInternodes);
    else
        % Gamma(shape k, scale theta): mean = k*theta, CV = 1/sqrt(k).
        k     = 1 / CV^2;
        theta = meanLength / k;

        lengths = theta * gammaRand(k, nInternodes);

        % Optional bounds via per-entry resampling (guarded against infinite loop).
        if minLength > 0 || isfinite(maxLength)
            for ii = 1:nInternodes
                attempts = 0;
                while lengths(ii) < minLength || lengths(ii) > maxLength
                    lengths(ii) = theta * gammaRand(k, 1);
                    attempts = attempts + 1;
                    if attempts > maxAttempts
                        error('sampleInternodeLengths:bounds', ...
                            ['Could not draw a length inside [%g, %g] after %d attempts. ', ...
                             'The bounds are likely in a near-zero-probability region for ', ...
                             'mean=%g, CV=%g.'], minLength, maxLength, maxAttempts, meanLength, CV);
                    end
                end
            end
        end
    end

    % ----- total-length clamp (also pins the mean exactly) -----
    if ~isempty(totalLength)
        lengths = lengths * (totalLength / sum(lengths));
        if (minLength > 0 || isfinite(maxLength)) ...
                && (min(lengths) < minLength || max(lengths) > maxLength)
            warning('sampleInternodeLengths:rescaleBounds', ...
                ['Rescaling to TotalLength pushed some lengths outside [%g, %g]. ', ...
                 'The total-length clamp takes precedence over the bounds.'], ...
                 minLength, maxLength);
        end
    end

    lengths = reshape(lengths, 1, []);   % ensure row vector
end


% ------------------------------------------------------------------------
function g = gammaRand(k, n)
%GAMMARAND  Gamma(shape=k, scale=1) samples using base MATLAB only.
%   Marsaglia & Tsang (2000) rejection method (uses RANDN and RAND, so no
%   Statistics Toolbox is required). Valid for any k > 0; for k < 1 it uses
%   the boosting identity  Gamma(k) = Gamma(k+1) * U^(1/k).
    k0    = k;
    boost = false;
    if k < 1
        boost = true;
        k = k + 1;
    end
    d = k - 1/3;
    c = 1 / sqrt(9 * d);

    g = zeros(1, n);
    i = 0;
    while i < n
        x = randn;
        v = (1 + c * x)^3;
        if v <= 0
            continue
        end
        u = rand;
        if log(u) < 0.5 * x^2 + d - d * v + d * log(v)
            i = i + 1;
            g(i) = d * v;
        end
    end

    if boost
        g = g .* (rand(1, n) .^ (1 / k0));
    end
end


% ------------------------------------------------------------------------
function val = getOption(args, name, default)
%GETOPTION  Simple name/value lookup (mirrors buildClampedAxon's helper).
    val = default;
    for kk = 1:2:numel(args) - 1
        if strcmpi(args{kk}, name)
            val = args{kk + 1};
        end
    end
end
