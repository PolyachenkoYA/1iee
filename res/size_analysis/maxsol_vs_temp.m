clear;
close all;

jobs = {'job111', 'job112', 'job113', 'job114'};
%jobs = {'job111', 'job113'};
%jobs = {'job111'};
stab_times = [200, 400, 600, 800, 1000];
%stab_times = [1000];
t = [1, 10, 15, 20, 25, 30, 35, 40];
draw_jobs = 0;
draw_stab_times = 1;

N_jobs = length(jobs);
N_st = length(stab_times);
N_lines = N_jobs * N_st;

if(draw_jobs)
    for j_i = 1:N_jobs
        fig = getFig('$T (C^\circ)$', 'maxsol/$N_{supercells}$', ['maxsol($T$ : $P$ = 1); ' jobs{j_i}]);
        for st_i = 1:N_st
            line_id = (j_i - 1) * N_st + (st_i - 1);
            line_name = [jobs{j_i} '_' num2str(stab_times(st_i))];
            draw_line(fig.ax, line_name, st_i, t);
        end
    end
end

if(draw_stab_times)
    for st_i = 1:N_st
        fig = getFig('$T (C^\circ)$', 'maxsol/$N_{supercells}$', ['maxsol($T$ : $P$ = 1); ' num2str(stab_times(st_i))]);
        for j_i = 1:N_jobs    
            line_id = (j_i - 1) * N_st + (st_i - 1);
            line_name = [jobs{j_i} '_' num2str(stab_times(st_i))];
            draw_line(fig.ax, line_name, j_i, t);
        end
    end
end


function draw_line(ax, line_name, line_id, t)
    maxsol_0 = load([line_name '.txt']);

    fit2 = polyfit(t, maxsol_0, 2);
    fit1 = polyfit(t, maxsol_0, 1);
    x_draw = linspace(min(t), max(t), 100);

    plot(ax, t, maxsol_0, 'o', ...
        'DisplayName', line_name, 'Color', getMyColor(line_id), 'LineWidth', 2);
    plot(ax, x_draw, polyval(fit1, x_draw), ...
        'HandleVisibility', 'off', 'Color', getMyColor(line_id), 'LineWidth', 1);
    %plot(ax, x_draw, polyval(fit2, x_draw), '--', ...
    %    'HandleVisibility', 'off', 'Color', getMyColor(line_id * 2));
end

