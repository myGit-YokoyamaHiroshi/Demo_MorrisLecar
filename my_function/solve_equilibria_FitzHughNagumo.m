function [v_eq, w_eq] = solve_equilibria_FitzHughNagumo(varargin)
    if length(varargin)==1    
        par  = varargin{1};
    else
        par  = varargin;
    end

    Iext = par{1};
    tau  = par{2};
    a    = par{3};
    b    = par{4};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the value of membrane potential at the equilibrium point
    %
    % Solve polynomial equations of the following form:
    %    coef(1) * v^3 + coef(2) * v^2 + coef(3) * v^1 + coef(4) * v^0 = 0
    coef = [1, 0, 1/b-1, -a/b - Iext];
    v_eq = roots(coef);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the value of recovery variable w at the equilibrium point
    w_eq = v_eq - v_eq.^3 + Iext;
    
    % exclude the complex value
    idx = find(imag(v_eq)==0 & imag(w_eq)==0 );

    v_eq = v_eq(idx);
    w_eq = w_eq(idx);
end