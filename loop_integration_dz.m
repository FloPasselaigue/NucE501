function int = loop_integration_dz(theta)
    % Compute the integrale of theta

    % Get the parameters needed in this particular workspace
    global s nLeg1 nLeg2 nCool nCells ;
    
    % derivative of z = f(s)
    dz_ds = zeros(nCells, 1);
    dz_ds(1:nLeg1) = 1/nLeg1 ;
    dz_ds(nLeg1+nCool+1:nLeg1+nCool+nLeg2) = -1/nLeg2 ;
    
    % integration matrix
    I = integration_matrix(s) ;
    
    int = sum( I * theta.*dz_ds );
    
end

function I = integration_matrix(s) 
    n = length(s) ;
    ds =s(2)-s(1) ;

    I = zeros(n,n) ;
    I(1,end) = 1 ;
    I(end,end) = 4 ;
    for i=1:n-1
        I(i,i) = 4 ;
        I(i+1,i) = 1 ;
        I(i,i+1) = 1 ;
    end
       
    I = I .* ds / 6 ;
end

