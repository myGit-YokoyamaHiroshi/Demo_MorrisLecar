function dXdt = MorrisLecar(X, varargin)
    V    = X(1);
    N    = X(2);
    
    if length(varargin)==1    
        par  = varargin{1};
    else
        par  = varargin;
    end
    
    C    = par{1};
    gL   = par{2};
    gK   = par{3};
    gCa  = par{4};
    VL   = par{5};
    VK   = par{6};
    VCa  = par{7};
    V1   = par{8};
    V2   = par{9};
    V3   = par{10};
    V4   = par{11};
    Iext = par{12};
    phi  = par{13}; 

    Minf = Sigm(V, V1, V2);
    Ninf = Sigm(V, V3, V4);
    
    dVdt = 1/C * (- gL  * (V - VL) ...
                  - gCa * Minf * (V - VCa) ...
                  - gK  * N  * (V - VK) + Iext);
    dNdt =  Lambda(V, V3, V4, phi) * (Ninf - N);

    dXdt = [dVdt, dNdt];
end

function val = Sigm(V, V1, V2)
    %%%% sigmoid function
    val =  1 / (1 + exp(-2 * (V - V1)/V2));
    % This function can be also expressed as: val = 0.5 * (1 + tanh((V - V1)/V2)); 
end

function lambda = Lambda(V, V1, V2, phi)
    lambda = phi * cosh((V-V1)/(2*V2));
end