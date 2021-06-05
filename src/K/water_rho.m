clear; close all;

T0 = 35;
T0 = 18;

TK2C = 273.15;
T = T0 + TK2C;
rho_sat = get_rho_sat(get_opt_rho_sat(TK2C, 0), T);

do_rho = 1;
do_L = 0;

w_extra0 = 1014 - 180;   % for T=35C
w_extra = [-40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]' + w_extra0;
w_extra = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610]' + w_extra0;
w_extra = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]' + w_extra0;
w_extra = [(0:50:450)'; w_extra];
%w_extra = (0:50:450)';

rho = [6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]';
rho = [6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8]';
rho = [5.76, 5.91, 6.77, 6.66, 7.38, 6.65, 6.4, 7.06, 7.17, 7.78]';
rho = [rho', 6.5, 9, 8, 6.1, 7, 10, 8, 7.85, 7, 8.5, 7.9, 8, 6.96, 6.36, 5.81, 7.27, 6.00, 7.70, 7.0, 8.0, 7.0, 7.86, 10, 7.3, 7.33, 8.0, 7.66, 8.43]';

d_rho = [0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]';
d_rho = [0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2]';
d_rho = [0.06, 0.06, 0.07, 0.065, 0.07, 0.07, 0.06, 0.07, 0.07, 0.07]';
%d_rho = [0.6, 1, 1, 1.2, 0.8, 0.7, 0.56, 1, 2, 1.5, 0.6, 1, 2, 0.3, 1, 3, 2, 0.5, 1, 1.5, 0.6, 2, 0.456, 0.348, 0.144, 0.800, 1.0, 1.0, 0.1, 1.5, 1.5, 1.18, 3.5, 1.5, 0.942, 0.4, 1.0, 0.182]';
d_rho = ones(size(rho)) * 0.4;

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
Nw = (w_extra + 207 * 8) * 2;
Na = 6.02 * 1e23;
Nprot = 8 * 2;
Vvac = 7.1 * 7.1 * 25 * 1e-27;
Vprot = 7.1 * 7.1 * (32.5 - 25) * 1e-27;
h = (Nw - rho * Vvac / 18 * Na) * 18 / (14530 * Nprot);
d_h = d_rho *  Vvac * Na / (14530 * Nprot);
v2 = 1 ./ (1 + h * 1.38);
d_v2 = 1.38 * d_h ./ (1 + h * 1.38).^2;
v1 = 1 - v2;
%x = Nw / Nprot * (1 / v1 - 1);
%V1 = Vprot / (Nw + x * Nprot);   % = Vprot / Nw * v1
V1 = 0.018 / (Na * 1000);   % rho_w = 1 g/cm^3
phi = rho / rho_sat;
d_phi = d_rho / rho_sat;
Kfit_x = (v2 .^ (-1/3) - 1) .* (5/3 * v2 .^ (-1/3) - 1) ./ v2.^2;
d_Kfit_x = abs((2*(20-28*v2.^(1/3)+9*v2.^(2/3)))./(9*v2.^(11/3))) .* d_v2;
Kfit_y = (log(phi ./ (1 - v2)) - v2) ./ v2.^2;
d_Kfit_y = sqrt((((((-2+v2).*v2)./(-1+v2)-2*log(phi./(1-v2)))./v2.^3) .* d_v2) .^ 2 + (d_phi ./ phi ./ v2.^2).^2);
Kfit_lin = polyfit(Kfit_x, Kfit_y, 1);
Kfit_K = Kfit_lin(1) / (V1 / (2 * T * 8.31 / Na));

fit_obj = fit(w_extra, rho, fittype('poly1'), 'Weight', 1 ./ d_rho);
linfit = coeffvalues(fit_obj);
d_linfit = confint(fit_obj);
d_linfit = (d_linfit(2, :) - d_linfit(1, :)) / 2;
d_rho_fnc = @(w)sqrt(d_linfit(2) ^ 2 + (d_linfit(1) * w) .^ 2);
w_ext_fnc = @(rho)([(rho - linfit(2)) / linfit(1),...
                    sqrt((rho * d_linfit(1) / linfit(1)^2).^2 + (d_linfit(2) / linfit(1))^2 + (d_linfit(1) * linfit(2) / linfit(1)^2)^2)]);

if(do_rho)
    fig_rho_Kfit = getFig('$(v_2^{-1/3} - 1) \cdot (5/3 \cdot v_2^{-1/3} - 1) / v_2^2$', '$(ln(P / P_0) - ln(1 - v_2) - v_2) / v_2^2$', '$y(x)$');
    errorbar(fig_rho_Kfit.ax, Kfit_x, Kfit_y, d_Kfit_y, d_Kfit_y, d_Kfit_x, d_Kfit_x, 'o', 'DisplayName', ['T = ' num2str(T0) '$C^{\circ}$']);
    plot(fig_rho_Kfit.ax, Kfit_x, polyval(Kfit_lin, Kfit_x), 'DisplayName', ['K = ' num2str(Kfit_K * 1e-6) ' MPa'])

    fig_rho_w = getFig('$w_{extra}$ ($H_2 O$ / cell)', '$\rho_w$ (g / $m^3$)', '$\rho_w(w_{extra})$');
    errorbar(fig_rho_w.ax, h, rho, d_rho, 'o', 'DisplayName', 'data');
    plot(fig_rho_w.ax, h, ones(size(w_extra)) * rho_sat, 'DisplayName', '$\rho_{sat}$');

    fig_rho_h = getFig('$\rho_w / \rho_{sat}$', '$h$ ($m_{w} / m_{prot}$)', '$h(\varphi)$');
    errorbar(fig_rho_h.ax, phi, h, d_h, d_h, d_phi, d_phi, 'o', 'DisplayName', ['T = ' num2str(T0) '$C^{\circ}$']);
end

if(do_L)
    fig_Lxy_w = getFig('$w_{extra}$ ($H_2 O$ / cell)', '$L$ (nm)', '$L(w_{extra})$');
    errorbar(fig_Lxy_w.ax, w_extra, Lx, d_Lx, 'o', 'DisplayName', '$L_x$');
    errorbar(fig_Lxy_w.ax, w_extra, Ly, d_Ly, 'o', 'DisplayName', '$L_y$');

    fig_Lxy_rho = getFig('$\rho_{vapor}$ (g/$m^3$)', '$L$ (nm)', '$L(w_{extra})$');
    errorbar(fig_Lxy_rho.ax, rho, Lx, d_Lx, d_Lx, d_rho, d_rho, 'o', 'DisplayName', '$L_x$');
    errorbar(fig_Lxy_rho.ax, rho, Ly, d_Ly, d_Ly, d_rho, d_rho, 'o', 'DisplayName', '$L_y$');
end
                
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


