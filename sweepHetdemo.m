% =========================================================================
% DEMO: FAST SMOKE TEST FOR STEP 2 HETEROGENEITY PIPELINE
% =========================================================================
clear; clc; close all;

fprintf('--- Starting Step 2 Heterogeneity Smoke Test ---\n');

% 1. Parameters for a FAST test
meanL     = 80;         % Mean internode length (um)
N_nodes   = 20;         % Small number of nodes for a fast solve
totLength = meanL * N_nodes; 

CV_levels = [0, 0.4, 0.8]; % Test Homogeneous (0), Mild (0.4), and High (0.8) variance
nReal     = 3;             % Only 3 realizations per level for speed

% 2. Sampler Bounds (CRUCIAL for stability)
% We cap the shortest internode at 15 um to prevent Crank-Nicolson ringing,
% and the max at 300 um to prevent massive truncation error matrices.
min_um = 15;  
max_um = 300; 

% 3. Run the Sweep
% We pass 'MaxSegments', 50 just for this demo so the solver doesn't 
% accidentally build a massive matrix if a rare 300um internode is drawn.
t_start = tic;

R_demo = sweepHetCVEnergy(meanL, CV_levels, ...
    'nReal', nReal, ...
    'TotalLength_um', totLength, ...
    'SamplerMin', min_um, ...
    'SamplerMax', max_um, ...
    'MaxSegments', 50, ...         % Cap resolution for speed
    'WindowFrac', [0.2 0.8], ...   % Measure interior 60%
    'Verbose', true, ...
    'Plot', true);                 % Auto-generate the plot

t_end = toc(t_start);
fprintf('\n--- Smoke Test Complete in %.1f seconds! ---\n', t_end);

% 4. Verify the Homogeneous Baseline (Acceptance Test)
% If CV_L = 0 works, the spread of CVs across realizations should be exactly 0.
cv0_idx = find(R_demo.CV_L == 0);
if ~isempty(cv0_idx)
    cv0_realizations = R_demo.cv(cv0_idx, :);
    disp('Acceptance Test: Homogeneous (CV=0) Velocities (should be identical):');
    disp(cv0_realizations);
end