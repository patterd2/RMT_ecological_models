tic;
%Main params
rng(1);
N=50;
mu = -1;
rho = -0.9;
%NB: The r parameters are a particular choice of a random set of response
%curves, and will be somewhat sensitive to the parameters.
r_mean = 0.5;
r_sd = 0.2;
r_spread = 2;

p_mean = 0;
p_sd = 0.5;
survival_threshold = 0.001;
alpha = 5;

i=1;
delta=1;


%Constructed params
rand_vec = randn(N,1);
%r = @(m)abs(r_mean + r_sd .* randn(N,1));
%Make r=1 at x=0.75, and r=-1 at x=0 and 1.5.
r = @(m)2*sech(r_spread*(m-(r_mean+r_sd.* rand_vec))).^2-1;
p = (p_mean + p_sd .* randn(1,N));
A = random_elliptic(N, mu, alpha, rho);
x0 = 0.5*ones(N,1);
m0 = 1;


[T, x, m] = Simulate_MetalLV(r, A, i, delta, p, x0, m0, 1e2);


close all;

% Create a figure
figure;

% Plot 1: Abundances
subplot(1, 3, 1);  % 1 row, 3 columns, position 1
plot(T,x,'linewidth',2);
title('Abundances');
xlabel('$t$','interpreter','latex');

% Plot 2: Metal Availability
subplot(1, 3, 2);  % position 2
plot(T,m,'linewidth',2);
title('Metal Availability');
xlabel('$t$','interpreter','latex');

% Plot 3: Growth rates
subplot(1, 3, 3);  % position 3
Ms = linspace(0,i*2,1e3);
plot(Ms, r(Ms), 'LineWidth',2);
title('Growth rates');
xlabel('$m$','interpreter','latex');


finalAbundance = x(end,:);
S_hat = sum(finalAbundance > survival_threshold);
prop_survived = S_hat / N

toc;

