% Fig. 3(a): uniqueness condition (9) holds -> two runs with different initial
% adoption converge to the same equilibrium.

rng(1);
n = 30;    % agents
m = 5;     % technologies
T = 350;   % horizon

%% Network and budget parameters (Assumption 1(i), 1(iii))
[W, tilde_W, ~, ~] = build_params_assumption(n, m);
beta  = rand(n,m);   beta  = beta  ./ max(1, sum(beta,2));
gamma = rand(n,m);   gamma = gamma ./ max(1, sum(gamma,2));

%% lambda, xi, delta chosen so that (9) holds
% The delta factor in (9) is >= 2 and grows like 1/delta^2, so delta is kept
% close to 1 to allow a non-negligible xi.
lambda = 0.30 * rand(n,m);
delta  = 0.45 + 0.3 * rand(n,m);              
xi     = 0.1 * rand(n,m) .* (1 - lambda);     % lambda + xi < 1

val = hyp_unique_value(lambda, xi, delta);
fprintf('[Case A] condition (9) LHS = %.4f  (< 1 required)\n', val);
assert(val < 1, 'condition (9) NOT satisfied (LHS = %.4f >= 1)', val);

%% Predisposition p = X0 > 0 (Assumption 1(ii)), same in both runs
X0 = rand(n,m);

%% Two different initial adoption conditions
A0_run1 = zeros(n,m);
A0_run1(:,1) = 0.1*rand(n,1);
A0_run1(:,2) = 0.05*rand(n,1);
A0_run1(:,5) = 0.1*rand(n,1);

A0_run2 = zeros(n,m);
A0_run2(:,3) = 0.15*rand(n,1);
A0_run2(:,4) = 0.1*rand(n,1);
A0_run2(:,1) = 0.05*rand(n,1);
A0_run2(:,5) = 0.3*rand(n,1);

%% Simulations
[~, A1, D1, X1] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0_run1, X0);
[~, A2, ~,  X2] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0_run2, X0);

meanA1 = squeeze(mean(A1,2));
meanA2 = squeeze(mean(A2,2));
meanD1 = squeeze(mean(D1,2));
meanX1 = squeeze(mean(X1,2));
meanX2 = squeeze(mean(X2,2));

fprintf('[Case A] Final adoption fractions (run 1): '); disp(meanA1(end,:));
fprintf('[Case A] Final adoption fractions (run 2): '); disp(meanA2(end,:));
fprintf('[Case A] Max adoption equilibrium difference: %.2e\n', max(abs(meanA1(end,:) - meanA2(end,:))));
fprintf('[Case A] Max opinion equilibrium difference : %.2e\n', max(abs(meanX1(end,:) - meanX2(end,:))));

%% Plot: one color per technology, run 1 solid, run 2 dashed
colors = lines(m);
figure; hold on;
h_run1 = zeros(1,m); h_run2 = zeros(1,m);
for k = 1:m
    h_run1(k) = plot(0:T-1, meanA1(:,k), '-',  'LineWidth', 2, 'Color', colors(k,:));
    h_run2(k) = plot(0:T-1, meanA2(:,k), '--', 'LineWidth', 2, 'Color', colors(k,:));
end

% Legend fills by rows: interleave handles so run 1 / run 2 form the two columns
leg1 = arrayfun(@(k) ['Adopters [', num2str(k), '], run 1'], 1:m, 'UniformOutput', false);
leg2 = arrayfun(@(k) ['Adopters [', num2str(k), '], run 2'], 1:m, 'UniformOutput', false);
handles = reshape([h_run1; h_run2], 1, []);
labels  = reshape([leg1;  leg2], 1, []);
legend(handles, labels, 'Location', 'best', 'NumColumns', 2, 'Interpreter','latex');
xlabel('Time', 'Interpreter','latex'); ylabel('Fraction of population','Interpreter','latex');
grid on; set(gca,'FontSize',16);
