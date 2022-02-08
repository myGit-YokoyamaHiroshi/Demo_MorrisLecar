function X_next = runge_kutta(X_now, dt, func, varargin)
    k1     = func(X_now, varargin);
    
    X_k2   = X_now + (dt/2) * k1;
    k2     = func(X_k2, varargin);
    
    X_k3   = X_now + (dt/2) * k2;
    k3     = func(X_k3, varargin);
    
    X_k4   = X_now + dt * k3;
    k4     = func(X_k4, varargin);

    X_next = X_now + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end