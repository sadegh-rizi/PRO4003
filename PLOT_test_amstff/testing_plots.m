%% Plot testing

% import data
cvFile = 'CVEnergy_table.csv';
precisionFile = 'results_aim1/precision.csv';

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
precision_sigma_t_ms = Tp.sigma_t_ms(:);
precision_sigma_per_mm = Tp.sigma_per_mm(:);


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

