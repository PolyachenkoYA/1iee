close all; clear;

%my_D = [1.248, 1.489, 1.7, 1.947, 2.228, 2.486, 2.81, 3.13, 3.45, 3.81, 4.18, 4.56];
%d_my_D = [0.003, 0.025, 0.02, 0.001, 0.025, 0.008, 0.01, 0.028, 0.003, 0.05, 0.04, 0.04];
my_t = [0.1, 5:5:55];

my_D = [2.077, 3.867, 8.346];
d_my_D = [0.004, 0.005, 0.02];
article_t = [273, 300, 350] - 273;
tip4p = [2.14, 3.75, 7.39];
d_tip4p = [0.04, 0.05, 0.1];
tip2005 = [1.24, 2.39, 5.37];
d_tip2005 = [0.03, 0.04, 0.09];

err = exp(abs(log(my_D ./ tip4p))) - 1;
err = my_D ./ tip4p - 1;
d_err = exp(abs(log(my_D ./ tip4p))) .* sqrt((d_my_D ./ my_D).^2 + (d_tip4p ./ tip4p).^2);
d_err = abs(my_D ./ tip4p) .* sqrt((d_my_D ./ my_D).^2 + (d_tip4p ./ tip4p).^2);

fig_D = getFig('$T$ ($C^{\circ}$)', '$D$ ($nm^2/ns$)', '$D(T)$');
errorbar(article_t, my_D, d_my_D, 'o', 'DisplayName', 'my');
errorbar(article_t, tip4p, d_tip4p, 'o', 'DisplayName', 'tip4p');
errorbar(article_t, tip2005, d_tip2005, 'o', 'DisplayName', 'tip4p/2005');

fig_err = getFig('$T$ ($C^{\circ}$)', '$D_{my} / D_{tip4p} -1$', '$err(T)$');
errorbar(article_t, err, d_err, 'o', 'DisplayName', 'my vs tip4p');

