function [M, N, P, dM_dw] = build_matrices(w_n, w_n_1)
    % Compute the matrices used in the linear step

    % Get the parameters needed in this particular workspace
    global nLeg1 nCool nHeat nCells V_Vh phi St_m ds dtau
    
    a_n_1 = w_n_1 * phi / (4 * ds) ;
    a_n   = w_n   * phi / (4 * ds) ;
    
    % Construction of M and N
    M = zeros(nCells, nCells);
    N = zeros(nCells, nCells);
    for i=1:nCells-1
        M(i,i) = 1/dtau ;
        M(i,i+1) = a_n_1 ;
        M(i+1,i) = -a_n_1 ;
        
        N(i,i) = 1/dtau ;
        N(i,i+1) = -a_n ;
        N(i+1,i) = a_n ;
    end
    
    % Cooler region
    M(nLeg1:nLeg1+nCool-1 , nLeg1:nLeg1+nCool-1) =... 
        M(nLeg1:nLeg1+nCool-1 , nLeg1:nLeg1+nCool-1) + St_m*eye(nCool) ;
    
    N(nLeg1:nLeg1+nCool-1 , nLeg1:nLeg1+nCool-1) =... 
        N(nLeg1:nLeg1+nCool-1 , nLeg1:nLeg1+nCool-1) - St_m*eye(nCool) ;
    
    M(end,end) = 1/dtau ;
    N(end,end) = 1/dtau ;
    
    
    % Construction of P
    P = zeros(nCells, nCells);
    P(end-nHeat+1 : end, end-nHeat+1 : end) = V_Vh * eye(nHeat) ;

    % Peridic boundary conditions
    M(1,end) = M(2,1);
    N(1,end) = N(2,1);
    P(1,end) = P(2,1);

    M(end,1) = M(1,2);
    N(end,1) = N(1,2);
    P(end,1) = P(1,2);
    
    % Construction of dM_dw
    dM_dw = M / w_n_1 ;
    for i=1:nCells
        dM_dw(i,i) = dM_dw(i,i)*w_n_1 ;
    end
    
end
    



























