% SpeedEnergyDemo - conduction velocity and Na+ energy cost of a cortical axon.
%
% Demonstrates the conduction-speed and energy modules on the
% Arancibia-Carcamo et al. (2017) cortical axon. Shows three things:
%   1. The one-call measurement pipeline (measureAxon).
%   2. The individual modules (conductionSpeed, energyConsumption).
%   3. A cross-check of conduction velocity against the shipped velocities()
%      helper.
%
% Run from anywhere:
%   >> SpeedEnergyDemo

% --- put the repository (and its subfolders) on the path ---
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
addpath(genpath(repoRoot));

% --- build the cortical axon parameter set ---
par = Carcamo2017CortexAxon();

%% 1. One-call pipeline ----------------------------------------------------
fprintf('Running cortical axon with ionic-current recording...\n');
m = measureAxon(par, 'Verbose', true, 'Vcross', -20);

fprintf('\n==== measureAxon: Carcamo2017 cortex ====\n');
fprintf('Internode length L         : %.1f um\n',        m.internodeLength_um);
fprintf('Conduction velocity        : %.3f m/s\n',       m.cv_m_per_s);
fprintf('Na+ charge per AP per mm   : %.4g nC/mm\n',     m.energy_nC_per_AP_per_mm);
fprintf('ATP per AP per mm          : %.4g ATP/mm\n',    m.atp_per_AP_per_mm);
fprintf('AP propagated across axon  : %d\n',             m.propagated);

%% 2. Individual modules, with full diagnostics ---------------------------
[V, IL, t, CURRENTS] = Model(par, [], false, 1);
dt = t(2) - t(1);

cvReg = conductionSpeed(V, IL, dt, 'Vcross', -20, 'Method', 'regression');
cv2pt = conductionSpeed(V, IL, dt, 'Vcross', -20, 'Method', 'twopoint');
en    = energyConsumption(CURRENTS, IL, 'NodeSegPerNode', par.geo.nnodeseg);

fprintf('\n==== conductionSpeed ====\n');
fprintf('Regression CV : %.3f m/s   (R^2 = %.5f over %d interior nodes)\n', ...
        cvReg.cv, cvReg.r2, numel(cvReg.nodesUsed));
fprintf('Two-point CV  : %.3f m/s\n', cv2pt.cv);

fprintf('\n==== energyConsumption ====\n');
fprintf('Na+ channels summed      : %s\n', strjoin(en.channelsUsed, ', '));
fprintf('Interior span            : %.3f mm over %d nodes\n', ...
        en.interiorLength_mm, numel(en.nodesUsed));
fprintf('Total Na+ charge per AP  : %.4g nC\n',     en.totalNaChargePerAP);
fprintf('Na+ charge per AP per mm : %.4g nC/mm\n',  en.chargePerAPperMM);
fprintf('Na+ ions  per AP per mm  : %.4g ions/mm\n', en.naIonsPerAPperMM);

%% 3. Cross-check CV against the shipped two-point velocities() helper -----
velMid = velocities(V, IL, dt, [20, 30], 'voltagecross', -20);
fprintf('\n==== cross-check ====\n');
fprintf('velocities() nodes 20->30: %.3f m/s\n', velMid);
fprintf('(should be close to the conductionSpeed estimates above)\n');
