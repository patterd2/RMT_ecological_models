function [T, x, m] = Simulate_MetalLV(r, A, i, delta, p, x0, m0, TEnd, dist, dist_start, dist_end, delta_i, enforce_death, survival_threshold)
%Simulate the MetalLV, with a slight tweak to prevent numerical blowup to
%negative infinity.
%U(1:N) are the N spevies, and U(end) is the metal.
N = length(x0);
if(dist)
    if(enforce_death)
        RHS = @(t,U)[U(1:N).*(r(U(end))-U(1:N)+A*U(1:N)).*(U(1:N)>survival_threshold)-U(1:N).*(U(1:N)<survival_threshold); i+delta_i*(t>dist_start)*(t<dist_end)-U(end).*(delta+p*U(1:N))];
    else
        RHS = @(t,U)[U(1:N).*(r(U(end))-U(1:N)+A*U(1:N)); i+delta_i*(t>dist_start)*(t<dist_end)-U(end).*(delta+p*U(1:N))];
    end
else
    RHS = @(t,U)[U(1:N).*(r(U(end))-U(1:N)+A*U(1:N)); i-U(end).*(delta+p*U(1:N))];
end
opts = odeset('reltol',1e-10,'abstol',1e-10,'maxstep',0.1,'NonNegative',1:N+1);
[T, U] = ode15s(RHS, linspace(0,TEnd,1e5),[x0;m0],opts);
x = U(:, 1:N); m = U(:,end);
end