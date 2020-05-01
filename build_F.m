function [F, dF_dx] = build_F(Re, theta_n, w_n, w_n_1)
    % Build the function F used in the non-linear step

    % Get the parameters needed in this particular workspace
    global dtau Gr p b D_L ;

    % functions f and g terms
    alpha = (p / (D_L * 4 * Re^b)) ;

    f     = @(x) x / dtau + alpha * x.^(2-b) ;  
    df_dx = @(x) 1 / dtau + (2 - b) * alpha * x.^(1-b) ;

    g     = @(x) x / dtau - alpha * x.^(2-b) ;
    dg_dx = @(x) 1 / dtau - (2 - b) * alpha * x.^(1-b) ;

    % integrated profiles terms
    [M, N, P, dM_dw] = build_matrices(w_n, w_n_1);
    one = ones(length(M(:,1)), 1) ;
    M_inv = M \ eye(length(M)) ;
    theta_n_1 = M_inv * (N*theta_n + P*one) ;
    dtheta_dw = - M_inv * dM_dw * theta_n_1 ;
    
    int_theta_n   = loop_integration_dz(theta_n);
    int_theta_n_1 = loop_integration_dz(theta_n_1);
    int_dtheta_dw = loop_integration_dz(dtheta_dw);
    
    F     = @(x) f(x) - g(w_n) - .5 * (Gr / Re^3) * (int_theta_n + int_theta_n_1) ;
    dF_dx = @(x)    df_dx(x)   - .5 * (Gr / Re^3) * int_dtheta_dw ;
end

















