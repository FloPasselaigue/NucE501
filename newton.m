function [x, iter] = newton(x0, f, df_dx, prec, it_max)
    % Use Newton's method comput find a zero of function f
    % starting from x0, with precision prec and a maximum  
    % number of iterations it_max
    
    x_old = x0 ;
    x_new = x_old - f(x_old) / df_dx(x_old) ;
    
    iter = 0 ;
    while (abs((x_old-x_new)/x_new) > prec) & (iter < it_max)
        x_old = x_new ;
        x_new = x_old - f(x_old) / df_dx(x_old) ;
        iter = iter + 1;
    end

    if(iter == it_max)
        warning('max iterations reached in Newton scheme; Residual: %s', num2str(f(x_new)));
    end
    
    x = x_new ;
end

