tic;
%Main params
rng(1);
N=10;
mu = -1;
rho = -0.9;
r_mean = 1;
r_sd = 0.5;
survival_threshold = 0.001;
alpha = 5;

%Constructed params
r = abs(r_mean + r_sd .* randn(N,1));
A = random_elliptic(N, mu, alpha, rho);
x0 = 0.5*ones(N,1);


[T, x] = Simulate_GLV(r, A, x0, 2e3);

plot(T,x);
axis tight;


finalAbundance = x(end,:);
S_hat = sum(finalAbundance > survival_threshold);
prop_survived = S_hat / N

toc;

