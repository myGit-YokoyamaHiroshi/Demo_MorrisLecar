function [v, w] = get_nullcline_FizHughNagumo(varargin)
    Iext = varargin{1};
    tau  = varargin{2};
    a    = varargin{3};
    b    = varargin{4};

    vmin = varargin{5};
    vmax = varargin{6};

    v    = linspace(vmin,vmax,100);

    w1   = v - v.^3 + Iext;
    w2   = 1/b .* (v - a);

    w    = [w1,; w2].';
end