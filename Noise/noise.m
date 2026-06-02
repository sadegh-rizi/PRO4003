%% Noise (Gausian noise)
clear; clc;

thisDirectory = fileparts(mfilename('fullpath'));
modelRoot = fileparts(thisDirectory);
addpath(genpath(modelRoot));

% diagnostics
% disp('Checking path:');
% which Model -all
% which Ford2015GBClat -all
% which GenerateEmptyParameterStructure -all
% which AddActiveChannel -all

%% Model and pars

% parameter structure
par0 = Carcamo2017OpticNerveAxon();
% Carcamo2017CortexAxon();
% Bakiri2011Cerebellum();
% Bakiri2011CorpusCallosum();
% Cullen2018CortexAxon();

% time
par0.sim.tmax.value = 10;               % simulation tmax
par0.sim.tmax.units = {1, 'ms', 1};     % units milisec
par0.sim.dt.value = 1;                  % 0.01-0.005 ms dt
par0.sim.dt.units = {1,'us',1};         % units microsec

% noise
par0.noise.amp.value = 0.05             % SD of noise current sample
par0.noise.amp.units = {1, 'nA', 1};    % nA
par0.noise.excludeStimNode = false;     % true = no noise at first stimulated node
par0.noise.seed = [123];                % [] = random each run


% Ls to test (change)
L_values = 30:10:170;

% variable allocation space (update)
velocities_ms = zeros(size(L_values));  % CVs
% metab_cost = ...

% Loop over Ls
for inter_l = 1:numel(L_values)
    L = L_values(inter_l);

    % build axon with intern L
    par = par0;                                 % start from baseline
    par = UpdateInternodeLength(par, L);        % homogns L = value
    %par = UpdateInternodeLength(par, Ls_vec);  % hetergns L = vector

    % run deterministic model (suppress verbose output)
    [Vm, internodeLength, tvec] = Model_noise(par, [], true);

    % measure and store simulation output
    dt    = tvec(2) - tvec(1);                  % time step in ms
    nodes = [20, 40];                           % measurement nodes
    % watch out... tricky
    nodes = min(nodes, size(Vm, 2));            % clamp to columns in Vm
    % --------------------
    vel   = velocities(Vm, internodeLength, dt, nodes);
    velocities_ms(inter_l) = vel;               % m/s

    fprintf('L = %d µm => velocity = %.2f m/s\n', L, vel);
end

%% Plotting
figure;
plot(L_values, velocities_ms, 'o-');
xlabel('Internode length (micromtr)');
ylabel('CV (m/s)');
title('Velocity vs internode L');