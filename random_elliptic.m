function Gamma = random_elliptic(N, mu, alpha, rho)
% Inputs:
%   N     - The dimension of the square matrix (N-by-N).
%   mu    - The mean of the distribution.
%   alpha - A scaling factor for the standard deviation.
%   rho   - correlation coefficient


% Coefficients for symmetric and antisymmetric parts.
c1 = sqrt((1 + rho) / 2);
c2 = sqrt((1 - rho) / 2);

X = randn(N);
Y = randn(N);

A = zeros(N);
for i = 1:N
    for j = (i+1):N
        % Assign A(i,j) and A(j,i) based on the elliptic law formula.
        x = X(i,j);
        y = Y(i,j);
        A(i,j) = c1 * x + c2 * y;
        A(j,i) = c1 * x - c2 * y;
    end
end

% Set the diagonal entries of A. These are independent standard
% normal variables.
A(1:N+1:end) = randn(N, 1);

Gamma = A / (alpha*sqrt(N)) + (mu / N);

end