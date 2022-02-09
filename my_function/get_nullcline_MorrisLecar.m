function [v, n] = get_nullcline_MorrisLecar(varargin)
    gL   = varargin{1};
    gK   = varargin{2};
    gCa  = varargin{3};
    VL   = varargin{4};
    VK   = varargin{5};
    VCa  = varargin{6};
    V1   = varargin{7};
    V2   = varargin{8};
    V3   = varargin{9};
    V4   = varargin{10};
    Iext = varargin{11};

    vmin = varargin{12};
    vmax = varargin{13};
    
    v    = linspace(vmin,vmax,100);

    Minf = Sigm(v, V1, V2);
    Ninf = Sigm(v, V3, V4);
    
    v_nullcline = (- gCa.*Minf.*(v-VCa) - gL.*(v-VL) + Iext)./(gK.*(v-VK));
    n_nullcline = Ninf;

    n = [v_nullcline,; n_nullcline].';
end


function val = Sigm(V, V1, V2)
    %%%% sigmoid function
    val =  1 ./ (1 + exp(-2 .* (V - V1)./V2));
    % This function can be also expressed as: val = 0.5 * (1 + tanh((V - V1)/V2)); 
end
