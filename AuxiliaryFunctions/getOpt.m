function val = getOpt(args, name, default)
%GETOPT  Standalone fallback for the getOpt() name/value helper.
%   Identical semantics to the local function of the same name at the bottom
%   of run_step2_production.m. This file exists only as a safety net: MATLAB
%   always prefers a function's own local/subfunction definition over a
%   path-level function of the same name, so this file is silently ignored
%   whenever run_step2_production.m's own local getOpt() resolves correctly.
%   It only gets called if something (a stale cached parse, a shadowing copy
%   of the file elsewhere on the path, etc.) stops that local resolution from
%   working -- in which case this restores the same behaviour instead of
%   throwing "Unrecognized function or variable 'getOpt'".
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name), val = args{k + 1}; end
end
end
