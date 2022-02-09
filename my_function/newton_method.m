function x = newton_method(x0, func, tol, max_iter, varargin)
    eps = 1E-10;
    
    for i = 1:max_iter
        func_diff = (func(x0 + eps, varargin) - func(x0 - eps, varargin))./(2*eps); 
        
        x = x0 - func(x0, varargin)./func_diff;

        if abs(x-x0)<=tol
            break
        end
        x0 = x;
    end
end