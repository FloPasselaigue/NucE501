function evolW = solveur(Re, initW, initT, withPlot)


    % Get the parameters needed in this particular workspace
    global nCells tau_max dtau ;   

        % Convergence criteria
    precNonLin = 1e-3 ;
    stepNonLin = 5 ;
    precNewton = 1e-4 ;
    iterMaxNewton = 30;
    erreur = false ;

        % Initial guess: we use the initial condition
    w_n = initW ;
    w_n_1 = w_n ;
    theta_n = initT ; 

        % Time loop initialization
    evolW = zeros(1, tau_max / dtau) ;
    evolW(1) = w_n ;
    iter = 1 ;

    idFlow = 2 ;
    if withPlot 
        prepFig2(idFlow, nCells, tau_max, dtau);
        h = plot(evolW, 'LineWidth', 2);
        pause(0);
    end

        % Time loop
    while iter * dtau < tau_max && ~erreur

        % ---------------- Non linear step: determine w_(n+1) ------------ 
        w_n_1_old = w_n_1 ;

        % Build the function
        [F, dF_dx] = build_F(Re, theta_n, w_n, w_n_1) ;

        % Find the root
        [w_n_1, ~] = newton(w_n_1, F, dF_dx, precNewton, iterMaxNewton);

        % Repeat until convergence
        iterNonLin = 1 ;
        while abs((w_n_1_old - w_n_1)/w_n_1_old) > precNonLin
            w_n_1_old = w_n_1 ;
            [F, dF_dx] = build_F(Re, theta_n, w_n, w_n_1) ;
            [w_n_1, ~] = newton(w_n_1, F, dF_dx, precNewton, iterMaxNewton);
           
            iterNonLin = iterNonLin + 1 ;
            if iterNonLin > stepNonLin 
                erreur = true ;
                break ;
            end
        end

        % -------------- Linear step: determine theta_(n+1) --------------

        % Build the matrices used to solve the temperature equation
        [M, N, P, ~] = build_matrices(w_n, w_n_1);
        one = ones(length(M(:,1)), 1) ;

        % Resolution
        theta_n_1 = M \ (N*theta_n + P*one) ;

        % ------------------------ Incrementation ------------------------
        iter = iter + 1 ;
        theta_n = theta_n_1 ;
        w_n = w_n_1 ;

        % --------------------- Result Plot/Storage  ---------------------
        evolW(iter) = w_n ;
   
        if withPlot && rem(iter,fix(iter/30)) == 0
            figure(idFlow)
            delete(h);
            h = plot(evolW, 'Color', 'k', 'LineWidth', 2);
            pause(0);
        end        
    end
    % end time loop
    
        % Final plot
    if withPlot
        figure(idFlow)
        delete(h);
        h = plot(evolW, 'Color', 'k', 'LineWidth', 2);
        pause(0);
    end
    
        % Error management
    if erreur
        evolW(end) = NaN ;
    end
        
    
end

% Auxiliary plot functions

function prepFig2(id, nCells, tau_max, dtau)
    figure(id) 
    clf
    set(gca, 'FontSize',18, 'YGrid', 'on')
    xlabel('Step');
    ylabel('Adimensionnal Mass Flow');
    title("Evolution of the mass flow - " +... 
        num2str(nCells) + " cells - \tau_m_a_x = " + ...
        num2str(tau_max) + " ; d\tau = " + num2str(dtau) );
    hold on
    pause(1);
end
