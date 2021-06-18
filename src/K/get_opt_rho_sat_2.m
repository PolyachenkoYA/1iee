function rho_sat = get_opt_rho_sat_2(h, rho, d_rho, verbose)
    %a0 = [max(rho) * 1.1, 1/3]';
    %a0 = [2.3, 1.7]';
    a0 = [3.5];
    a_opt = fminunc(@(a)get_err(a, h, rho, d_rho), a0);
    rho_sat = a_opt(1);    
    
    if(verbose)
        getFig('$\rho / \rho_{sat}$', '$h$');
        errorbar(rho / rho_sat, h, [], [], d_rho / rho_sat, d_rho / rho_sat, 'DisplayName', 'data');
        
        x_draw = linspace(min(rho), max(rho), 100);
        plot(x_draw / rho_sat, get_h_th(a_opt, x_draw), 'DisplayName', ['$\rho_0 = ' num2str(rho_sat) '$']);
        plot(x_draw / rho_sat, get_h_th(a0, x_draw), 'DisplayName', 'initial fit');        
    end
    
end

function h = get_h_th(a, rho)
    %rho0 = a(1);
    %c = a(2) - 1;
    %phi = rho / rho0;
    %h = phi .* ((c - 1) ./ (1 + (c - 1) * phi) + 1 ./ (1 - phi));
    h = 1 ./ (1 - rho / a(1)) - 1;
end

function err = get_err(a, h, rho, d_rho)
    err = sum(((h - get_h_th(a, rho))).^2);
end