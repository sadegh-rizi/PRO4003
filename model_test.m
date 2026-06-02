
% Load parameter set
par = Cullen2018CortexAxon(); % needs to be replaced with baseline model


% noise
par.noise.amp.value = 0.05;             % SD of noise current sample
par.noise.amp.units = {1, 'nA', 1};    % nA
%par.noise.excludeStimNode = false;     % true = no noise at first stimulated node
%par.noise.seed = 123;                % [] = random each run


% internode lengths
internode_lengths = [40, 60, 100]; % internode lengths that needs to be tested

% Run model
[mean_delay, std_delay, failure_rate, norm_jitter, delays] = precision_sim(par, internode_lengths, 2);
