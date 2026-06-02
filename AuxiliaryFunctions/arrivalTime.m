function tc = arrivalTime(Vtrace, dt, vcross, doInterp)
%ARRIVALTIME  Time of the first upward crossing of a voltage threshold (ms).
%
%   tc = ARRIVALTIME(Vtrace, dt, vcross, doInterp)
%
%   Returns the time (ms) at which a node's membrane potential first rises
%   through vcross, with optional sub-sample linear interpolation; NaN if the
%   trace never crosses (propagation failure / no spike).
%
%   This is the single shared definition of "arrival time" used by both
%   conductionSpeed() (deterministic conduction velocity) and the precision
%   sweep (per-trial distal arrival time under noise), so the two readouts
%   can never drift apart.
%
%   Inputs:
%       Vtrace   - (T x 1) membrane potential of one node over time (mV).
%       dt       - sample interval (ms). Sample n is at time (n-1)*dt.
%       vcross   - threshold-crossing voltage (mV), e.g. -20.
%       doInterp - true (default) for sub-sample linear interpolation of the
%                  crossing time; false to return the nearest sample time.

VariableDefault('doInterp', true);

Vcol  = Vtrace(:);
above = Vcol > vcross;
if ~any(above)
    tc = NaN;
    return
end

m = find(above, 1, 'first');     % first sample at/above threshold
if m == 1 || ~doInterp
    tc = (m - 1) * dt;           % sample n has time (n-1)*dt
    return
end

% Linear interpolation between samples m-1 (below) and m (at/above).
v0 = Vcol(m - 1);
v1 = Vcol(m);
denom = v1 - v0;
if denom == 0
    frac = 0;
else
    frac = (vcross - v0) / denom;
end
frac = min(max(frac, 0), 1);
tc = ((m - 1) - 1) * dt + frac * dt;
end
