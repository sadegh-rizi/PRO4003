function val = getOption(args, name, default)
%GETOPTION  Standalone fallback for the getOption() helper used locally in
%   sweepHetCVEnergy.m and sweepHetPrecision.m (see getOpt.m for why this
%   file exists -- it's a safety net, not the primary definition).
val = default;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name), val = args{k + 1}; end
end
end
