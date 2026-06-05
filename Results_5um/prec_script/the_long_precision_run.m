%% Precission sweep ======================¿===============================
% 100 trials
% 0.02 nA noise
% dt 0.25 us

intrnode_ls = 30:5:160;

R2 = sweepPrecision(intrnode_ls, 'NoiseAmp_nA', 0.02, 'nTrials', 100 , ...
    'BaseSeef', 1000, 'Tmax', 8, 'Dt_us', 0.25, 'SaveCsv', 'precision.csv');

% ========================================================================