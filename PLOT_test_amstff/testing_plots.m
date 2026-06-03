%% Plot testing

% import data
cvFile = 'CVEnergy_table.csv';
precisionFile = 'results_aim1/precision_batch_120-140-160.csv';

Tcv = readtable(cvFile);
Tp = readtable(precisionFile);

% extract vectors from Energy_CV
L_um = Tcv.L_um(:);
L_actual_um = Tcv.L_actual_um(:);           % actual length after clamped build, um
CV_mps = Tcv.cv_m_per_s(:);                 % conduction velocity, m/s
Na_pC_mm = Tcv.charge_pC_per_mm(:);         % Na+ charge proxy, pC/AP/mm
ATP_mol_mm = Tcv.atp_mol_per_mm(:);         % mol ATP/AP/mm
propagated = logical(Tcv.propagated(:));

% extract vectors from precision (dummy values)
precision_sigma_t_ms = [ ...
    0.095293207, ...
    0.096545629, ...
    0.097798050, ...
    0.099050471, ...
    0.100302893, ...
    0.101555314, ...
    0.102807736, ...
    0.104060157, ...
    0.105312579, ...
    0.106565000, ... % real, L = 120 um
    0.126088000, ...
    0.145611000, ... % real, L = 140 um
    0.135946000, ...
    0.126281000  ... % real, L = 160 um
    ];
precision_sigma_per_mm = [ ...
    0.028915997, ... 30
    0.029312442, ... 40
    0.030708887, ... 50
    0.029105331, ... 60
    0.030501776, ... 70
    0.031098221, ... 80
    0.030294666, ... 90
    0.031691110, ... 100
    0.036087555, ... 110
    0.032484000, ... % real, L = 120 um
    0.038612500, ... 130
    0.044741000, ... % real, L = 140 um
    0.040987500, ... 150
    0.037234000  ... % real, L = 160 um
    ];


%% 2D scatter plot CV vs Energy cost for diff internode Ls

figure('Name','2D speed-energy tradeoff','Color','w');
scatter(CV_mps, ATP_mol_mm, 80, L_um, 'filled');
hold on;
plot(CV_mps, ATP_mol_mm, 'k-', 'LineWidth', 0.75);  % links sweep order
for i = 1:numel(L_um)
    text(CV_mps(i), ATP_mol_mm(i), sprintf('  %g', L_um(i)), 'FontSize', 8);
end
hold off;

colormap(jet(256));
cb = colorbar;
ylabel(cb, 'Internode length L (\mum)');
xlabel('Conduction velocity (m/s)');
ylabel('ATP cost (mol ATP / AP / mm)');
title('CV vs ATP cost by internode length');
grid on; box on;


%% 3D mesh plot (experimental and prolly a bit broken)
figure('Name','3D CV-energy-precision tradeoff','Color','w');
hold on;

X = CV_mps;
Y = ATP_mol_mm;
Z = precision_sigma_per_mm;
C = L_um;

% Build mesh
try
    tri = delaunay(X, Y);
    trisurf(tri, X, Y, Z, C, ...
        'FaceAlpha', 0.35, ...
        'EdgeColor', [0.35 0.35 0.35]);
catch meshErr
    warning('Could not build trisurf mesh: %s. Showing curve/scatter only.', meshErr.message);
end

% colored points by intrnode L
scatter3(X, Y, Z, 70, C, 'filled');

% Label points with L
for i = 1:numel(C)
    text(X(i), Y(i), Z(i), sprintf('  %g', C(i)), 'FontSize', 8);
end

hold off;

% --
colormap(jet(256));
cb = colorbar;
ylabel(cb, 'Internode length L (\mum)');

xlabel('Conduction velocity (m/s)');
ylabel('ATP cost (amol ATP / AP / mm)');
zlabel('Timing precision \sigma / distance (ms/mm)');

title('Interactive CV-ATP-precision tradeoff; color = internode length');
grid on; box on;
view(135, 30);
rotate3d on;

