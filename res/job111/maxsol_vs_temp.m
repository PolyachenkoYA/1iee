clear;
close all;

t = [274, 283, 288, 293, 298, 303, 308, 313] - 273;
maxsol_0 = [1127, 1092.75, 1081.2, 1071.44, 1049.3, 1035.8, 1023, 1002.8];
t_pred = [1, 10, 15, 20, 25, 30, 35, 40];

fit2 = polyfit(t, maxsol_0, 2);
fit1 = polyfit(t, maxsol_0, 1);
x_draw = linspace(min(t), max(t), 100);
maxsol_pred2 = polyval(fit2, t_pred);
maxsol_pred1 = polyval(fit1, t_pred);

fig = getFig('$T (C^o)$', 'maxsol/$N_{supercells}$', 'maxsol(T)');
plot(fig.ax, t, maxsol_0, 'o', 'DisplayName', '1x1x1', 'Color', getMyColor(1));
plot(x_draw, polyval(fit2, x_draw), 'DisplayName', 'linear fit', 'Color', getMyColor(1));
%plot(x_draw, polyval(fit2, x_draw), '--', 'DisplayName', '$x^2$ fit', 'Color', getMyColor(2));
disp(maxsol_pred1 - 1000);

