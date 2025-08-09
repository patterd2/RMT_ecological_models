function [T, x, m] = Simulate_MetalLV(r, A, i, delta, p, x0, m0, TEnd)
    %Simulate the MetalLV, with a slight tweak to prevent numerical blowup to
    %negative infinity.
    %U(1:N) are the N spevies, and U(end) is the metal.
    N = length(x0);
    RHS = @(t,U)[U(1:N).*(r(U(end))-U(1:N)+A*U(1:N)); i-U(end).*(delta+p*U(1:N))];
    opts = odeset('reltol',1e-10,'abstol',1e-10,'maxstep',0.1,'NonNegative',1:N+1);
    [T, U] = ode15s(RHS, linspace(0,TEnd,1e4),[x0;m0],opts);
    x = U(:, 1:N); m = U(:,end);
end