function lengths = samplePoisson(nInternodes, meanLength, minLength, maxLength)
%SAMPLEPOISSON Generate heterogeneous internode lengths using a Poisson process.
%   LENGTHS = SAMPLEPOISSON(NINTERNODES, MEANLENGTH) returns a row vector of
%   NINTERNODES positive internode lengths drawn from an exponential
%   distribution with mean MEANLENGTH..
%
%   LENGTHS = SAMPLEPOISSON(NINTERNODES, MEANLENGTH, MINLENGTH, MAXLENGTH)
%   allows you to specify optional minimum and maximum allowable
%   internode lengths.  Internodes whose length falls outside the
%   [MINLENGTH, MAXLENGTH] range are resampled until they satisfy the
%   constraints. 
%
%   IMPORTANT; internode number must be precomputed after axon length 
%   is clamped, otherwise use the nInternodes argument to manually define 
%   the number of internodes in the axon.

%   Example usage:
%     % Generate 50 internodes with mean length 60 µm and default bounds
%     lengths = samplePoisson(50, 60);
%     % Plot histogram of generated lengths
%     histogram(lengths);

%   See also RANDOM, EXPRND.

    narginchk(2,4);
    % Set default bounds if they are not provided
    if nargin < 3 || isempty(minLength)
        minLength = 1; % µm – smallest plausible node of Ranvier
    end
    if nargin < 4 || isempty(maxLength)
        maxLength = inf; % no upper limit by default
    end
    if ~isscalar(nInternodes) || nInternodes <= 0 || floor(nInternodes) ~= nInternodes
        error('nInternodes must be a positive integer');
    end
    if meanLength <= 0
        error('meanLength must be a positive scalar');
    end
    if minLength <= 0
        error('minLength must be positive');
    end
    if maxLength < minLength
        error('maxLength must be greater than or equal to minLength');
    end
    % Preallocate output
    lengths = zeros(1, nInternodes);
    % Generate lengths one by one, resampling as needed
    for ii = 1:nInternodes
        len = 0;
        while len < minLength || len > maxLength
            % Exponential distribution with mean equal to meanLength
            % Equivalent to inter-event intervals of a Poisson process
            u = rand();
            % Avoid log(0)
            if u == 0
                u = eps;
            end
            len = -meanLength * log(u);
        end
        lengths(ii) = len;
    end
end
