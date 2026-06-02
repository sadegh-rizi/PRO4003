
% Load parameter set
par = Cullen2018CortexAxon(); % needs to be replaced with baseline model

% internode lengths
internode_lengths = []; % internode lengths that needs to be tested

% Run model
[mean_delay, std_delay, failure_rate, norm_jitter, delays] = precision_sim(par, internode_lengths, 5);
