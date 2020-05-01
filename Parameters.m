function Parameters()
    %% Reactor geometry https://doi.org/10.13182/NT10-74
    global l1 lCool l2 lHeat H V_Vh phi D_L D_h_core A_core;
    lHeat = 2.5 ;
    lCool = 2.1 ;
    H = 4.8 ;
    l1 = H + .5 * (lCool - lHeat) ;
    l2 = H + .5 * (lHeat - lCool) ;
    
    D_core = .95 ;
    D_pin = 14e-3 ;
    H_pin = 5.2 ;
    P_D = 1.1  ;
    nAssemblies = 18 ;
    nPinAssembly = 169 ;

    % hydraulic diameter of interior channel
    A_int = D_pin^2 * ((sqrt(3)/4) * P_D^2 - pi/8) ;
    P_int = pi * D_pin / 2 ;
    D_h_int = 4 * A_int / P_int ;

    % hydraulic diameter of core
    A_core = nAssemblies * nPinAssembly * 2 * A_int ; % .5 pin / channel
    P_core = nAssemblies * nPinAssembly * P_int + pi*D_core ;
    D_h_core = 4 * A_core / P_core ;

    D_L = D_h_core / lHeat ;
    V_Vh = 8 ; % 
    phi = 3 ; % 

    
    %% Discretization
    global ds s nLeg1 nLeg2 nCool nHeat nCells tau_max dtau ;
    
    % Number of cells
    nLeg1 = floor(1 / ds) ;
    nLeg2 = floor(1 / ds) ;
    nCool = floor(lCool / ds / H) ;
    nHeat = floor(lHeat / ds / H) ;
    
    % Space vector
    nCells = nLeg1 + nLeg2 + nCool + nHeat ;
    lTot = 2 * H + lCool + lHeat ;
    S = linspace(0, lTot, nCells) ;
    s = S / H ; % adimensionnal  
    
    % Time 
    tau_max = 30 ;
    dtau    = .005;
    
    %% Coolant properties (liquid Sodium @430K) 
    % https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    global Cp mu k beta rho Pr;
    
    Cp = 1.22e3 ;
    mu = 6.54e-4 ;
    k = 89.36 ;
    beta = 2.44e-4 ; 
    rho = 912 ;

    Pr = Cp * mu / k ;
    
    %% Adimensionnal numbers
    global St_m Gr p b ;
    p    = 22.26 ;
    b    = 1 ;
    Gr   = 1e10 ;
    St_m = 15 ; 
    
    
end






















