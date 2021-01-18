close all; clear;

my_t = [0.1, 5:5:55];

%my_D = [2.077, 3.867, 8.346];
%d_my_D = [0.004, 0.005, 0.02];
my_D = [2.0833, 2.3548, 2.6759, 2.9911, 3.3278, 3.7376, 4.1991, 4.4784, 4.8772, 5.3572, 5.7396, 6.21449];
d_my_D = [0.03, 0.03, 0.016, 0.015, 0.04, 0.03, 0.04, 0.03, 0.01,  0.08, 0.016,  0.11];
D_dens05 = [0.35090293463501687, 0.3729193788028815, 0.4277494751544869, 0.4567176592438946, 0.5444280443490108, 0.5360602169627409, 0.6196245598786887, 0.729197075860696, 0.6800428508710554, 0.7792474018336464, 0.7924840824005273, 0.8382789124881586];
D_cum09 = [0.28374138899914786, 0.33315631653812156, 0.3815744748952225, 0.43994058315972034, 0.5082903766817644, 0.5998387954725222, 0.7104885895349945, 0.7909502414405744, 0.8473407341031498, 0.9746638667280892, 1.0417477166480824, 1.1875705978025932];

article_t = [273, 300, 350] - 273;
tip4p = [2.14, 3.75, 7.39];
d_tip4p = [0.04, 0.05, 0.1];
tip2005 = [1.24, 2.39, 5.37];
d_tip2005 = [0.03, 0.04, 0.09];

%err = exp(abs(log(my_D ./ tip4p))) - 1;
%err = my_D ./ tip4p - 1;
%d_err = exp(abs(log(my_D ./ tip4p))) .* sqrt((d_my_D ./ my_D).^2 + (d_tip4p ./ tip4p).^2);
%d_err = abs(my_D ./ tip4p) .* sqrt((d_my_D ./ my_D).^2 + (d_tip4p ./ tip4p).^2);

fig_D = getFig('$T$ ($C^{\circ}$)', '$D$ ($nm^2/ns$)', '$D(T)$');
errorbar(my_t, my_D, d_my_D, 'o', 'DisplayName', 'my GROM tip4p');
plot(my_t, D_dens05, 'o', 'DisplayName', '$p_{dens} \approx 0.5$ $ns/nm^2$');
plot(my_t, D_cum09, 'o', 'DisplayName', '$P_{cum} \approx 0.9$');
errorbar(article_t, tip4p, d_tip4p, 'o', 'DisplayName', 'tip4p');
errorbar(article_t, tip2005, d_tip2005, 'o', 'DisplayName', 'tip4p/2005');

fig_cmp = getFig('$T$ ($C^{\circ}$)', '$y$', '$D_{prot}$ vs $D_{bulk}$');
plot(my_t, D_cum09 ./ my_D, 'o', 'DisplayName', '$D_{cum} / D_{bulk}$');
plot(my_t, D_dens05 ./ my_D, 'o', 'DisplayName', '$D_{dens} / D_{bulk}$');

%fig_err = getFig('$T$ ($C^{\circ}$)', '$D_{my} / D_{tip4p} -1$', '$err(T)$');
%errorbar(article_t, err, d_err, 'o', 'DisplayName', 'my vs tip4p');

