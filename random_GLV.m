%% Random gLV with multiple simulations over alpha (interaction strength)
% Parameters
N = 100;
T_max = 20;
dt = 0.1;
seed = 123;
survival_threshold = 0.001;
L = 20;
enforce_nonneg = true;
% Model parameters
r_mean = 1;
r_sd = 0.2;
mu = 0;
xi = 0.0;
init_abund = 0.5;
% Alpha values to loop over
alpha_values = 1:0.1:3;
n_alpha = numel(alpha_values);
% Results storage
avg_prop_survived_vec = zeros(n_alpha,1);
avg_m_hat_vec = zeros(n_alpha,1);
avg_convergence_vec = zeros(n_alpha,1);              % Added for convergence
convergence_threshold = 1e-5;                        % Relative change for convergence

%tic
for a_idx = 1:n_alpha
    alpha = alpha_values(a_idx);
    all_props = zeros(L,1);
    all_mhats = zeros(L,1);
    all_convergences = zeros(L,1);                   % Added for convergence
    for ell = 1:L
        % Intrinsic growth rates
        r = abs(r_mean + r_sd .* randn(N,1)); % nonnegative growth rates
        % Random interaction matrix
        A = random_elliptic(N, mu, alpha, xi);
        % Initial condition
        if isscalar(init_abund)
            x0 = init_abund * ones(N,1);
        else
            if numel(init_abund) ~= N
                error('Initial abundance vector must have length N.');
            end
            x0 = init_abund(:);
        end

        % Simulate the model in time
        [T, X] = Simulate_GLV(r, A, x0, T_max);

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

        % Convergence check (based on relative change between last two time points)
        if size(X, 1) > 1
            prevAbundance = X(end-1, :);
            rel_change = abs(finalAbundance - prevAbundance) ./ max(abs(finalAbundance));
            converged_species = (rel_change < convergence_threshold);
            prop_converged = sum(converged_species)/N;
        else
            prop_converged = NaN;
        end
        all_convergences(ell) = prop_converged;
    end
    % Store average metrics for this alpha
    avg_prop_survived_vec(a_idx) = mean(all_props);
    avg_m_hat_vec(a_idx) = mean(all_mhats, 'omitnan');
    avg_convergence_vec(a_idx) = mean(all_convergences, 'omitnan');   % Save average convergence
end

%% Plotting
% plot settings
lw = 2.5;              % Line width
ms = 8;                % Marker size
fontSize = 20;         % Font size

% Calculate standard deviations
std_prop_survived_vec = std(all_props, 0, 1); % Per-alpha std (L simulations)
std_m_hat_vec = std(all_mhats, 0, 1, 'omitnan'); % Per-alpha std (L simulations, omit NaN)

% Errorbar plot: Proportion survived vs alpha
figure('Name','Survival vs Alpha');
errorbar(alpha_values, avg_prop_survived_vec, std_prop_survived_vec, 'o-', ...
    'LineWidth', lw, ...
    'MarkerSize', ms, ...
    'Color', [0.2 0.4 0.8], ...
    'MarkerEdgeColor', [0.2 0.4 0.8], ...
    'MarkerFaceColor', [0.2 0.4 0.8]);
xlabel('\alpha (interaction strength)', 'FontSize', fontSize);
ylabel('Average proportion survived', 'FontSize', fontSize);
title('Species survival vs \alpha (mean \pm std)', 'FontSize', fontSize+2);
grid on;
set(gca, 'FontSize', fontSize);
box on;
axis tight;

% Errorbar plot: m_hat vs alpha
figure('Name','m_hat vs Alpha');
errorbar(alpha_values, avg_m_hat_vec, std_m_hat_vec, 'o-', ...
    'LineWidth', lw, ...
    'MarkerSize', ms, ...
    'Color', [0.85 0.33 0.1], ...
    'MarkerEdgeColor', [0.85 0.33 0.1], ...
    'MarkerFaceColor', [0.85 0.33 0.1]);
xlabel('\alpha (interaction strength)', 'FontSize', fontSize);
ylabel('m (mean square of surviving species)', 'FontSize', fontSize);
title('Community abundance (m) vs \alpha (mean \pm std)', 'FontSize', fontSize+2);
grid on;
set(gca, 'FontSize', fontSize);
box on;
axis tight;
