function J = jacobian_matrix_FitzHughNagumo(X, varargin)
    v    = X(1);
    w    = X(2);

    if length(varargin)==1    
        par  = varargin{1};
    else
        par  = varargin;
    end

    Iext = par{1};
    tau  = par{2};
    a    = par{3};
    b    = par{4};

    J    = [1 - 3*v.^2,      -1;
                1./tau, -b./tau];
end