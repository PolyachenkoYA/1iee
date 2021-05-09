clear; close all;

w_extra = [-40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]';
w_extra = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610]';
w_extra = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]';

rho = [6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]';
rho = [6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8]';
rho = [6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8, 6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]';

d_rho = [0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]';
d_rho = [0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2]';
d_rho = [0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2, 0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]';

d_rho_min = 0.75;

wrong_ids = d_rho < d_rho_min;
wght = 1 ./ d_rho;

fit_obj = fit(w_extra, rho, fittype('poly1'), 'Weight', 1 ./ d_rho);
linfit = coeffvalues(fit_obj);
d_linfit = confint(fit_obj);
d_linfit = (d_linfit(2, :) - d_linfit(1, :)) / 2;
d_rho_fnc = @(w)sqrt(d_linfit(2) ^ 2 + (d_linfit(1) * w) .^ 2);
w_ext_fnc = @(rho)([(rho - linfit(2)) / linfit(1),...
                    sqrt((rho * d_linfit(1) / linfit(1)^2).^2 + (d_linfit(2) / linfit(1))^2 + (d_linfit(1) * linfit(2) / linfit(1)^2)^2)]);

fig = getFig('$w_{extra}$ ($H_2 O$ / chain)', '$\rho_w$ (g / $cm^3$)', '$\rho_w(w_{extra})$');
errorbar(fig.ax, w_extra, rho, d_rho, 'o', 'DisplayName', 'data');
plot(fig.ax, w_extra, polyval(linfit, w_extra), ...
    'DisplayName', 'linfit', 'Color', getMyColor(2));
plot(fig.ax, w_extra, polyval(linfit, w_extra) + d_rho_fnc(w_extra), '--', ...
    'HandleVisibility', 'off', 'Color', getMyColor(2));
plot(fig.ax, w_extra, polyval(linfit, w_extra) - d_rho_fnc(w_extra), '--', ...
    'HandleVisibility', 'off', 'Color', getMyColor(2));

% yt = w_ext_fnc(polyval(linfit, w_extra));
% plot(fig.ax, yt(:, 1), polyval(linfit, w_extra), ...
%     'DisplayName', 'test', 'Color', getMyColor(3));
% plot(fig.ax, yt(:, 1) - yt(:, 2), polyval(linfit, w_extra), ...
%     'HandleVisibility', 'off', 'Color', getMyColor(3));
% plot(fig.ax, yt(:, 1) + yt(:, 2), polyval(linfit, w_extra), ...
%     'HandleVisibility', 'off', 'Color', getMyColor(3));


