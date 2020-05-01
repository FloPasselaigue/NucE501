function delta = check_stability(evolW, n)
    % Compute the evolution of the standard deviation

    l = length(evolW) / n ;
    dev = zeros(1, l);
    
    for i=1:l
        dev(i) = std( evolW(1 + (i-1)*n : i*n ) );
    end

    delta = (dev(2) - dev(1)) / dev(1) ;
end