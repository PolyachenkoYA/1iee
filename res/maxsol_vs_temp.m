clear;
close all;

jobs = {'job111', 'job112', 'job113'};
jobs = {'job111', 'job113'};
stab_times = [200, 400, 600, 800, 1000];
stab_times = [400, 600];
t = [274, 283, 288, 293, 298, 303, 308, 313] - 273;

N_jobs = length(jobs);
N_st = length(stab_times);
N_lines = N_jobs * N_st;
%maxsol_0 = [[1127.0418749292355, 1091.7501653791155, 1081.2097465600436, 1071.4355133779354, 1049.3102655592338, 1035.793411911533, 1023.0620698732068, 1002.8045724860353];
%            [1111.3541417260426, 1091.0404142445884, 1074.50431774713, 1058.92041946324, 1041.6509775793527, 1029.6817636250005, 1010.713405803926, 1000.6829731284107];
%            [1111.61188756838, 1088.6194872415094, 1072.5419372408371, 1059.5612634716779, 1040.1654093279215, 1027.6061852774062, 1011.4730049685628, 994.9958354159417]];

t_pred = [1, 10, 15, 20, 25, 30, 35, 40];

fig = getFig('$T (C^\circ)$', 'maxsol/$N_{supercells}$', 'maxsol($T$ : $P$ = 1)');
for j_i = 1:N_jobs
    for st_i = 1:N_st
        line_id = (j_i - 1) * N_st + (st_i - 1);
        line_name = [jobs{j_i} '_' num2str(stab_times(st_i))];
        maxsol_0 = load([line_name '.txt']);
        
        fit2 = polyfit(t, maxsol_0, 2);
        fit1 = polyfit(t, maxsol_0, 1);
        x_draw = linspace(min(t), max(t), 100);
        maxsol_pred2 = polyval(fit2, t_pred);
        maxsol_pred1 = polyval(fit1, t_pred);

        plot(fig.ax, t, maxsol_0, 'o', ...
            'DisplayName', line_name, 'Color', getMyColor(line_id), 'LineWidth', 2);
        plot(x_draw, polyval(fit1, x_draw), ...
            'HandleVisibility', 'off', 'Color', getMyColor(line_id), 'LineWidth', 1);
        %plot(x_draw, polyval(fit2, x_draw), '--', ...
        %    'HandleVisibility', 'off', 'Color', getMyColor(line_id * 2));
        %disp(maxsol_pred1 - 1000);
    end
end

