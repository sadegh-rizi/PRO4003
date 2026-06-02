function [mean_delay, std_delay, failure_rate, norm_jitter, delays] = precision_sim(baseline_par, internode_lengths, n_trials)
% This code runs 100 simulations to calculate delay, std, failure rate and
% normalized jitter per internode length.

% INPUT: baseline parameters in a .m file, batch of internode lengths to simulate
% OUTPUT: mean delay, std delay, failure rate and normalized jitter per
% internode length.

% (requires: set_internode_length and arrival_time)


n_lengths = length(internode_lengths);
n_nodes = baseline_par.geo.nnode;

% empty vectors to store results per internode length
mean_delay = nan(n_lengths,1);
std_delay = nan(n_lengths,1);
failure_rate = nan(n_lengths,1);
norm_jitter = nan(n_lengths,1);
delays = nan(n_trials,n_lengths);

%%% SIMULATION %%%
% loop internode length
for il = 1 : n_lengths

    L_i = internode_lengths(il);
    L_total_um    = L_i * n_nodes; % total internode length
    fprintf('Running internode length = %d µm ...\n', L_i);

    % set the internode length with pre-existing function
    par_modified = UpdateInternodeLength(baseline_par, L_i);

    % store delays of all trials for this internode length
    %delays = nan(n_trials,1);
    
    % run 100 simulations and store conduction delay (if AP occurred)
    for trial =  1: n_trials
        % run the model
        [MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = Model_noise(par_modified); 
        
        % determine conduction delay
        tarrival_first = arrival_time_detector(MEMBRANE_POTENTIAL(:,1), TIME_VECTOR);
        tarrival_last = arrival_time_detector(MEMBRANE_POTENTIAL(:,end), TIME_VECTOR);
        conduction_delay = tarrival_last - tarrival_first; % from first to last AP
        delays(trial,il) = conduction_delay;
    end
   
    % calculate the results for this internode length
    valid = ~isnan(delays(:,il));
    n_valid = sum(valid);

    mean_delay(il) = mean(delays(valid,il)); %% change for index
    std_delay(il) = std(delays(valid,il));
    norm_jitter(il) = std_delay(il) / L_total_um; % normalize the std by dividing by axon length

    failure_rate(il) = 1 - n_valid / n_trials;

    fprintf('  mean delay = %.3f ms | jitter = %.4f ms | norm jitter = %.4f ms | failures = %.0f%%\n', ...
            mean_delay(il), std_delay(il), norm_jitter(il), failure_rate(il)*100);
end
end