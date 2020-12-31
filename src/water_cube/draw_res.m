clear; close all;

data = load('times_data/4times_12models.dat')';
t = data(end, :);
K_data = data(1 : end-1, :);
N_models = size(K_data, 1);
K_flucts = mean(K_data);
d_K_flucts = std(K_data) / sqrt(N_models);

K_dV = [2.083, 2.41, 1.925, 2.009] * 1e9;
d_K_dV = [0.2, 0.412, 0.247, 0.5] * 1e9 / 1.5;

getFig('$t_{stab}$ (ns)', '$K$ (GPa)', 'water cube, $K(t_{stab})$, $T = 30$ $C^\circ$');
errorbar(t, K_dV * 1e-9, d_K_dV * 1e-9, 'o', 'DisplayName', 'fluct method', 'LineWidth', 1.5);
errorbar(t, K_flucts * 1e-9, d_K_flucts * 1e-9, 'o', 'DisplayName', 'dV method', 'LineWidth', 1.5);
plot([0, max(t)], [1, 1] * 2.1, '--', 'DisplayName', 'known', 'Color', [0,0,0], 'LineWidth', 2);
