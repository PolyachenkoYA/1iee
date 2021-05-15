clear; close all;

TK2C = 273.15;
T0 = 35;
rho_sat = get_rho_sat(get_opt_rho_sat(TK2C, 0), T0 + TK2C);

w_extra0 = 1014 - 180;   % for T=35C
w_extra = [-40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]' + w_extra0;
w_extra = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610]' + w_extra0;
w_extra = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]' + w_extra0;
w_extra = [(0:50:450)'; w_extra];
w_extra = (0:50:450)';

rho = [6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]';
rho = [6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8]';
rho = [5.23, 5.6, 5.9, 6.13, 6.43, 6.33, 6.0, 8.12, 9, 7.8]';
%rho = [5.04, 5.6, 5.5, 7, 6.3, 6.4, 6.25, 7.85, 8.2, 7, 6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8, 6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]';

d_rho = [0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]';
d_rho = [0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2]';
d_rho = [0.13, 1.6, 1.3, 1.5, 1.3, 1.4, 0.76, 2.04, 2.6, 2.1]';
%d_rho = [0.6, 1, 1, 1.2, 0.8, 0.7, 0.56, 1, 2, 1.5, 0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2, 0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]';

Ldata = [[7.3304770548, 0.015432846612721736, 7.1722916656, 0.018129885433487338]; ...
        [7.181661062399999, 0.018338095855549075, 7.1532357532, 0.01483626308301014]; ...
        [7.154355831600001, 0.024599868982025933, 7.1407179884, 0.019213137300676995]; ...
        [7.203556148800001, 0.01737954576049267, 7.157062537200001, 0.019548653034002526]; ...
        [7.14633481, 0.016616372016788138, 7.199240510800001, 0.01718012919565751]; ...
        [7.172980352399999, 0.01559303956660839, 7.1008510712, 0.018091184556565953]; ...
        [7.137558040800001, 0.03651943615975383, 7.1998307732, 0.018970557328981178]; ...
        [7.0711723064, 0.018759809651308277, 7.1621753956, 0.01655031349844167]; ...
        [7.072170213600001, 0.017807973590410975, 7.178682769999999, 0.020135832643014787]; ...
        [7.1689189096, 0.016244877394693623, 7.1224542972, 0.016652257547109712]];
Lx = Ldata(:, 1);
d_Lx = Ldata(:, 2);
Ly = Ldata(:, 3);
d_Ly = Ldata(:, 4);

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

fig_Lxy_w = getFig('$w_{extra}$ ($H_2 O$ / cell)', '$L$ (nm)', '$L(w_{extra})$');
errorbar(fig_Lxy_w.ax, w_extra, Lx, d_Lx, 'o', 'DisplayName', '$L_x$');
errorbar(fig_Lxy_w.ax, w_extra, Ly, d_Ly, 'o', 'DisplayName', '$L_y$');

fig_Lxy_rho = getFig('$\rho_{vapor}$ (g/$m^3$)', '$L$ (nm)', '$L(w_{extra})$');
errorbar(fig_Lxy_rho.ax, rho, Lx, d_Lx, d_Lx, d_rho, d_rho, 'o', 'DisplayName', '$L_x$');
errorbar(fig_Lxy_rho.ax, rho, Ly, d_Ly, d_Ly, d_rho, d_rho, 'o', 'DisplayName', '$L_y$');
                
fig_rho_w = getFig('$w_{extra}$ ($H_2 O$ / cell)', '$\rho_w$ (g / $m^3$)', '$\rho_w(w_{extra})$');
errorbar(fig_rho_w.ax, w_extra, rho, d_rho, 'o', 'DisplayName', 'data');
plot(fig_rho_w.ax, w_extra, ones(size(w_extra)) * rho_sat, 'DisplayName', '$\rho_{sat}$');
% plot(fig.ax, w_extra, polyval(linfit, w_extra), ...
%     'DisplayName', 'linfit', 'Color', getMyColor(2));
% plot(fig.ax, w_extra, polyval(linfit, w_extra) + d_rho_fnc(w_extra), '--', ...
%     'HandleVisibility', 'off', 'Color', getMyColor(2));
% plot(fig.ax, w_extra, polyval(linfit, w_extra) - d_rho_fnc(w_extra), '--', ...
%     'HandleVisibility', 'off', 'Color', getMyColor(2));

% yt = w_ext_fnc(polyval(linfit, w_extra));
% plot(fig.ax, yt(:, 1), polyval(linfit, w_extra), ...
%     'DisplayName', 'test', 'Color', getMyColor(3));
% plot(fig.ax, yt(:, 1) - yt(:, 2), polyval(linfit, w_extra), ...
%     'HandleVisibility', 'off', 'Color', getMyColor(3));
% plot(fig.ax, yt(:, 1) + yt(:, 2), polyval(linfit, w_extra), ...
%     'HandleVisibility', 'off', 'Color', getMyColor(3));


