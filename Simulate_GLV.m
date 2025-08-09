function [T, x] = Simulate_GLV(r, A, x0, TEnd)
    %Simulate the GLV, with a slight tweak to prevent numerical blowup to
    %negative infinity
    RHS = @(t,x)x.*(r-x+A*x);
    opts = odeset('RelTol',1e-10, 'AbsTol',1e-10,'NonNegative',1:length(x0));
    [T, x] = ode15s(RHS, linspace(0,TEnd,1e4),x0,opts);
end