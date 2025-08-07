% Matrix size
N = 1000;
% Bias term
mu = 0;
% variance scaling parameter
alpha = 1;
% Correlation parameter (-1 < rho < 1)
rho = 0.25;

Gamma = random_elliptic(N,mu,alpha,rho);

% Compute eigenvalues
eigvals = eig(Gamma);

% Find eigenvalue with largest real part
[~, idx_max_real] = max(real(eigvals));
eigval_max = eigvals(idx_max_real);
eigvals_main = eigvals;
eigvals_main(idx_max_real) = [];

% Plot eigenvalues
figure;
scatter(real(eigvals_main), imag(eigvals_main), 25, 'filled', 'MarkerFaceColor', [0.2 0.2 0.8]); hold on;
scatter(real(eigval_max), imag(eigval_max), 40, 'filled', 'MarkerFaceColor', 'r');

% Plot ellipse boundary
a = 1 + rho; 
b = 1 - rho; 

% Parametrize the ellipse
t = linspace(0, 2*pi, 100); 
x = a * cos(t)/alpha; 
y = b * sin(t)/alpha; 

% Plot the ellipse
plot(x, y,'--k','LineWidth',2);

xlabel('Real Part'); ylabel('Imaginary Part');
title(['Eigenvalues of \Gamma (bias \mu = ', num2str(mu),...
    ', correlation \rho = ', num2str(rho), ', \alpha = ', num2str(alpha),'), N = ', num2str(N)]);
legend('Eigenvalues', 'leading eigenvalue', 'Elliptic Law');
grid on;
