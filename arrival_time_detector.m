function t_last_node = arrival_time_detector(membrane_potential, time_vector)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
threshold = -20; % mV


below  = membrane_potential(1:end-1) <  threshold;
    above  = membrane_potential(2:end)   >= threshold;
    crossings = find(below & above);

    if isempty(crossings)
        t_last_node = NaN;   % AP failed to propagate 
    else
        % Linear interpolation for sub-timestep precision
        i  = crossings(1);
        dt = time_vector(i+1) - time_vector(i);
        dV = membrane_potential(i+1)  - membrane_potential(i);
        t_last_node = time_vector(i) + dt * (threshold - membrane_potential(i)) / dV;
    end
end
