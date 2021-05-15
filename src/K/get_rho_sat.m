function rho = get_rho_sat(a, T)
    %rho = exp(a(1) - a(2) ./ (T + a(3))) ./ T;
    rho = exp(a(1) - a(2) ./ T) ./ T.^a(3);
end
