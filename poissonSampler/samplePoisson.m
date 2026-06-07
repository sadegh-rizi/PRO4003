function lengths = samplePoisson(nInternodes, meanLength_um, minLength_um, maxLength_um, seed)
%SAMPLEPOISSON Sample heterogeneous internode lengths from a bounded Poisson distribution.
%
%   lengths = samplePoisson(nInternodes, meanLength_um, minLength_um, maxLength_um)
%   returns a 1 x nInternodes vector of internode lengths in micrometres.
%
%   The sampled values follow a Poisson distribution with mean meanLength_um,
%   truncated to the interval [minLength_um, maxLength_um].
%
%   Example:
%       L = samplePoisson(50, 80, 40, 160);
%       par = UpdateInternodeLength(par, L);

    narginchk(4, 5);

    if nargin == 5 && ~isempty(seed)
        rng(seed);
    end

    % input checks --------------------------------------------------------
    if ~isscalar(nInternodes) || nInternodes <= 0 || floor(nInternodes) ~= nInternodes
        error('nInternodes must be a positive integer.');
    end

    if ~isscalar(meanLength_um) || meanLength_um <= 0
        error('meanLength_um must be a positive scalar.');
    end

    if ~isscalar(minLength_um) || ~isscalar(maxLength_um)
        error('minLength_um and maxLength_um must be scalar values.');
    end

    if minLength_um <= 0 || maxLength_um <= 0
        error('minLength_um and maxLength_um must be positive.');
    end

    if maxLength_um < minLength_um
        error('maxLength_um must be greater than or equal to minLength_um.');
    end

    % Poisson values are integer-valued.
    minLength_um = ceil(minLength_um);
    maxLength_um = floor(maxLength_um);

    if maxLength_um < minLength_um
        error('No integer lengths exist inside the requested min/max range.');
    end

    % bounded Poisson probability mass function (PMF) ---------------------
    k = minLength_um:maxLength_um;

    % formula log(P(X = k)) = k*log(lambda) - lambda - log(k!)
    % gammaln(k+1) is log(k!) and is numerically stable.
    logP = k .* log(meanLength_um) - meanLength_um - gammaln(k + 1);

    % Shift before exponentiating to avoid numerical underflow.
    P = exp(logP - max(logP));

    % Truncate and renormalise.
    P = P / sum(P);

    CDF = cumsum(P);

    % sample from the cummulative distribution function (CDF) -------------
    lengths = zeros(1, nInternodes);

    for ii = 1:nInternodes
        u = rand();
        idx = find(u <= CDF, 1, 'first');
        lengths(ii) = k(idx);
    end
end