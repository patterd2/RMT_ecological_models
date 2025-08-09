%% Random generalized Lotka-Volterra with multiple simulations over alpha

% Parameters
N = 100;
T = 100;
dt = 0.1;
seed = 123;
survival_threshold = 0.001;
L = 20;
enforce_nonneg = true;

% Model parameters
r_mean = 1;
r_sd = 0.2;
mu = 0;
rho = 0.0;
init_abund = 0.5;

% Alpha values to loop over
alpha_values = 1:0.1:3.5;
n_alpha = numel(alpha_values);

% Results storage
avg_prop_survived_vec = zeros(n_alpha,1);
avg_m_hat_vec = zeros(n_alpha,1);

%tic

for a_idx = 1:n_alpha
    alpha = alpha_values(a_idx);

    all_props = zeros(L,1);
    all_mhats = zeros(L,1);

    for ell = 1:L
        rng(seed + ell);  % reproducibility with variation

        % Intrinsic growth rates
        r = abs(r_mean + r_sd .* randn(N,1)); % nonnegative growth rates

        % Random interaction matrix
        A = random_elliptic(N, mu, alpha, rho);

        % Initial condition
        if isscalar(init_abund)
            x0 = init_abund * ones(N,1);
        else
            if numel(init_abund) ~= N
                error('Initial abundance vector must have length N.');
            end
            x0 = init_abund(:);
        end


        % NOTE THIS IS TURNED OFF IN FAVOUR OF BETTER NUMERICS AT ABOUT 10X
        % SPEED COST; IF UNCOMMENTED, ALSO NEED TO SWAP THE DIMENSIONS OF X
        % AROUND TO MAKE THE finalAbundance CALCULATION CORRECT
        % % Time integration
        % tvec = 0:dt:T;
        % nt = numel(tvec);
        % X = zeros(N, nt);
        % X(:,1) = x0;
        %
        % for k = 2:nt
        %     x = X(:,k-1);
        %     dx = x .* r - x.^2 + x .* (A * x);
        %     x_new = x + dt * dx;
        %     if enforce_nonneg
        %         x_new(x_new < 0) = 0;
        %     end
        %     X(:,k) = x_new;
        % end

        % Simulate the model in time
        [T, X] = Simulate_GLV(r, A, x0, 100);

        % Metrics
        finalAbundance = X(end,:);
        S_hat = sum(finalAbundance > survival_threshold);
        prop_survived = S_hat / N;
        if S_hat > 0
            m_hat = sum(finalAbundance.^2) / S_hat;
        else
            m_hat = NaN;
        end

        all_props(ell) = prop_survived;
        all_mhats(ell) = m_hat;
    end

    % Store average metrics for this alpha
    avg_prop_survived_vec(a_idx) = mean(all_props);
    avg_m_hat_vec(a_idx) = mean(all_mhats, 'omitnan');

    %fprintf('alpha = %.2f: avg prop_survived = %.4f, avg m_hat = %.4f\n', ...
    %    alpha, avg_prop_survived_vec(a_idx), avg_m_hat_vec(a_idx));
end

%% Plotting

% plot settings
lw = 2.5;              % Line width
ms = 8;                % Marker size
fontSize = 16;         % Font size

% Plot: Proportion survived vs alpha
figure('Name','Survival vs Alpha');
plot(alpha_values, avg_prop_survived_vec, '-o', ...
    'LineWidth', lw, ...
    'MarkerSize', ms, ...
    'Color', [0.2 0.4 0.8], ...
    'MarkerEdgeColor', [0 0 0], ...
    'MarkerFaceColor', [0.2 0.4 0.8]);
xlabel('\alpha (interaction strength)', 'FontSize', fontSize);
ylabel('Average proportion survived', 'FontSize', fontSize);
title('Species survival vs \alpha', 'FontSize', fontSize+2);
grid on;
set(gca, 'FontSize', fontSize);
box on;
axis tight;

% Plot: m_hat vs alpha
figure('Name','m\_hat vs Alpha');
plot(alpha_values, avg_m_hat_vec, '-s', ...
    'LineWidth', lw, ...
    'MarkerSize', ms, ...
    'Color', [0.85 0.33 0.1], ...
    'MarkerEdgeColor', [0 0 0], ...
    'MarkerFaceColor', [0.85 0.33 0.1]);
xlabel('\alpha (interaction strength)', 'FontSize', fontSize);
ylabel('m (mean square of surviving species)', 'FontSize', fontSize);
title('Community abundance (m) vs \alpha', 'FontSize', fontSize+2);
grid on;
set(gca, 'FontSize', fontSize);
box on;
axis tight;

