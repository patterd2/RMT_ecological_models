%% Sweep alpha, plot proportion of species survived (colored by convergence status)

tic

T_max = 200;
N = 100;
survival_threshold = 0.001;
mu = 0.0;
xi = 0;
p_mean = 0.0;
p_sd = 0.0;
i = 1;
delta = 0.4;

alpha_values = 1:0.2:3;
n_alpha = numel(alpha_values);
n_sims = 25;
rel_threshold = 1e-3; % Relative change for convergence

% Preallocate arrays
all_alphas = repmat(alpha_values, n_sims, 1); % n_sims x n_alpha
all_props = zeros(n_sims, n_alpha);           % n_sims x n_alpha
all_mean_abund = zeros(n_sims, n_alpha);      % Mean final abundance
convergence_pass = false(n_sims, n_alpha);    % n_sims x n_alpha (logical)

for idx = 1:n_alpha
    alpha = alpha_values(idx);
    for sim = 1:n_sims
        % Generate random coefficients for r
        p1 = unifrnd(-0.15, -0.05, N, 1);
        p2 = unifrnd(0.2, 0.4, N, 1);
        p3 = unifrnd(0.2 , 0.9, N, 1);
        r = @(m) p1.*m.^2 + p2.*m + p3;
        p = (p_mean/N + (p_sd/sqrt(N)) .* randn(1,N));
        A = random_elliptic(N, mu, alpha, xi);

        % Random initial conditions
        x0 = unifrnd(0.5, 1.5, N, 1);
        m0 = unifrnd(1, 6, 1);

        [T, x, m] = Simulate_MetalLV(r, A, i, delta, p, x0, m0, T_max);

        finalAbundance = x(end, :);
        S_hat = sum(finalAbundance > survival_threshold);
        prop_survived = S_hat / N;
        all_props(sim, idx) = prop_survived;

        % Mean final abundance (all species)
        all_mean_abund(sim, idx) = sum(finalAbundance.^2)/S_hat;

        % Convergence check
        if size(x,1) > 1
            prevAbundance = x(end-1, :);
            rel_change = abs(finalAbundance - prevAbundance) ./ max(abs(finalAbundance), eps); % avoid zero division
            converged_species = (rel_change < rel_threshold);
            prop_converged = sum(converged_species)/N;
            convergence_pass(sim, idx) = prop_converged > 0.975; % e.g., declare whole run passed if >95% converged
        else
            convergence_pass(sim, idx) = false; % cannot check with only one step
        end
    end
end

% Flatten data for scatter plot
all_alphas_flat = all_alphas(:);
all_props_flat = all_props(:);
all_abundance_flat = all_mean_abund(:);
convergence_pass_flat = convergence_pass(:);

% Color setup: green for pass, red for fail
scatter_colors = repmat([0.6,0,0], numel(all_alphas_flat), 1); % default: dark red
scatter_colors(convergence_pass_flat,:) = repmat([0,0.4,0], sum(convergence_pass_flat), 1);

%% First figure: Survivor proportions
figure;
hold on;
scatter(all_alphas_flat, all_props_flat, 40, scatter_colors, 'o', 'LineWidth',1.25); % open circles, colored
xlabel('Interaction Strength \alpha','FontWeight','bold','FontSize',13);
ylabel('Proportion of Species Survived','FontWeight','bold','FontSize',13);
title('Survival Proportion','FontWeight','bold','FontSize',14);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
grid on;
box on;
hold off;

%% Second figure: Mean final abundance
figure;
hold on;
scatter(all_alphas_flat, all_abundance_flat, 40, scatter_colors, 'o', 'LineWidth',1.25); % open circles, colored
xlabel('Interaction Strength \alpha','FontWeight','bold','FontSize',13);
ylabel('Mean Final Abundance','FontWeight','bold','FontSize',13);
title('Mean Square of Final Abundance','FontWeight','bold','FontSize',14);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
grid on;
box on;
hold off;

%% Plot average proportion of survivors

% Calculate mean and std of survival proportion for each alpha
mean_surv = mean(all_props, 1);        % 1 x n_alpha
std_surv = std(all_props, 0, 1);       % 1 x n_alpha

% Plot mean survival proportion vs alpha, with error bars
figure;
errorbar(alpha_values, mean_surv, std_surv, 'o-', ...
    'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0.2 0.8], 'MarkerFaceColor', 'auto');
xlabel('Interaction Strength \alpha', 'FontWeight', 'bold', 'FontSize', 13);
ylabel('Average Proportion Survived', 'FontWeight', 'bold', 'FontSize', 13);
title('Mean \pm SD of Survival Proportion', 'FontWeight', 'bold', 'FontSize', 14);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
grid on;
box on;

%% Plot average of mean square abundances

% Calculate mean and std of mean square abundances for each alpha
mean_msq_abundance = mean(all_mean_abund, 1, 'omitnan');   % 1 x n_alpha
std_msq_abundance = std(all_mean_abund, 0, 1, 'omitnan');  % 1 x n_alpha

% Plot mean square final abundance vs alpha, with error bars
figure;
errorbar(alpha_values, mean_msq_abundance, std_msq_abundance, 's-', ...
    'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.8 0.2 0.1], 'MarkerFaceColor', 'auto');
xlabel('Interaction Strength \alpha', 'FontWeight', 'bold', 'FontSize', 13);
ylabel('Mean Square Final Abundance', 'FontWeight', 'bold', 'FontSize', 13);
title('Mean \pm SD of Final Square Abundance', 'FontWeight', 'bold', 'FontSize', 14);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
grid on;
box on;


toc;
