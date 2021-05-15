function a_opt = get_opt_rho_sat(TK2C, verbose)
    T = [300, 325, 350, 375]';
    rho = [5.806, 24.79, 83.05, 231.3]';
    d_rho = [1.531, 3.992, 14.42, 66.32]' / 100;

    cut_ind = 4;
    T(cut_ind:end) = [];
    rho(cut_ind:end) = [];
    d_rho(cut_ind:end) = [];

    T0 = TK2C + 35;

    a0 = [(log(rho(1)*T(1))*T(1) - log(rho(2)*T(2))*T(2)) / (T(1) - T(2)), ...
          (log(rho(1)*T(1)) - log(rho(2)*T(2))) / (1/T(2) - 1/T(1)), ...
          1]';
    a_opt = fminunc(@(a)get_err(a, T, rho, d_rho), a0);

    if(verbose)
        getFig('T (C)', '$\rho_{dat}$ (g/$m^3$)', ['$\rho_{sat}(T); \rho(' num2str(T0 - TK2C) ' C^{\circ}) = ' num2str(rho_sat(a_opt, T0)) '$ (g/$m^3$)']);
        errorbar(T - TK2C, rho, d_rho, 'o', ...
                 'DisplayName', 'data');
        T_draw = linspace(min(T), max(T), 100);
        plot(T_draw - TK2C, get_rho_sat(a_opt, T_draw), ...
             'DisplayName', 'fit');
        % plot(T_draw - TK2C, get_rho_sat(a0, T_draw), ...
        %      'DisplayName', 'initial');
    end
end

function err = get_err(a, T, rho, d_rho)
    err = sum(((get_rho_sat(a, T) - rho) ./ d_rho).^2, 'all');
end
