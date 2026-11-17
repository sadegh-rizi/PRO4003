function tf = hasOpt(args, name)
%HASOPT  Standalone fallback for the hasOpt() helper (see getOpt.m for why
%   this file exists -- it's a safety net, not the primary definition).
tf = false;
for k = 1:2:numel(args) - 1
    if strcmpi(args{k}, name), tf = true; return; end
end
end
