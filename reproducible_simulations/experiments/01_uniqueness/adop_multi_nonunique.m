%%%%%%%%%%%%%%%%%%%% ADOZIONE multi-tecnologia — Caso B
%% Regime in cui l'ipotesi di unicita' (eq:hyp-unique) NON e' soddisfatta
%% -> l'equilibrio finale dipende dalla condizione iniziale di adozione.
rng(3);
n = 30;    % numero di agenti
m = 5;     % numero di tecnologie
T = 600;   % orizzonte temporale

%% Parametri di rete e budget (Assumption (i), (iv))
[W, tilde_W, ~, ~] = build_params_assumption(n, m);
beta  = rand(n,m);        gamma  = 0.5*rand(n,m);
beta  = beta  ./ max(1, sum(beta,2));   % vicino al limite Assumption(iv): sum_k beta_i^k <= 1
gamma = gamma ./ max(1, sum(gamma,2));  % anello di retroazione di A forte

%% lambda, xi, delta calibrati per VIOLARE eq:hyp-unique
lambda = 0.20 + 0.35 * rand(n,m);                    % lambda_i^k in [0.60, 0.95) -> denom. piccolo
xi     = (0.90 + 0.09 * rand(n,m)) .* (1 - lambda);  % xi vicino al limite (1-lambda), lambda_i^k+xi_i^k<1
delta  = 0.02 + 0.06 * rand(n,m);                    % delta piccolo -> RHS grande

val = hyp_unique_value(lambda, xi, delta);
fprintf('[Caso B] eq:hyp-unique LHS = %.4f  (>= 1 atteso)\n', val);
assert(val >= 1, 'eq:hyp-unique risulta soddisfatta (LHS = %.4f < 1): questo non e'' il Caso B', val);

%% Opinioni iniziali strettamente positive (Assumption (iii)), fissate uguali nei due run
X0 = 0.1+0.8*rand(n,m);

%% Due diverse condizioni iniziali di ADOZIONE (non delle opinioni) -- stesse dei due run del Caso A
A0_run1 = zeros(n,m);
A0_run1(:,1) = 0.8*rand(n,1);
A0_run1(:,2) = 0.5*rand(n,1);
A0_run1(:,5) = 0.5*rand(n,1);

A0_run2 = zeros(n,m);
A0_run2(:,3) = 0.6*rand(n,1);
A0_run2(:,4) = 0.3*rand(n,1);
A0_run2(:,1) = 0.2*rand(n,1);

%% Simulazioni
[~, A1, ~, X1] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0_run1, X0);
[~, A2, ~, X2] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0_run2, X0);

meanA1 = squeeze(mean(A1,2));
meanA2 = squeeze(mean(A2,2));
meanX1 = squeeze(mean(X1,2));
meanX2 = squeeze(mean(X2,2));

fprintf('[Caso B] Frazioni finali adozione (run 1): '); disp(meanA1(end,:));
fprintf('[Caso B] Frazioni finali adozione (run 2): '); disp(meanA2(end,:));
fprintf('[Caso B] Differenza massima equilibrio adozione: %.2e\n', max(abs(meanA1(end,:) - meanA2(end,:))));
fprintf('[Caso B] Differenza massima equilibrio opinioni : %.2e\n', max(abs(meanX1(end,:) - meanX2(end,:))));

%% Plot confronto: un colore distinto per ogni tecnologia (stesso colore nelle
%% due run), run 1 = linea continua ( - ), run 2 = linea tratteggiata (--).
colors = lines(m);                              % colori distinti per tecnologia
figure; hold on;
h_run1 = zeros(1,m); h_run2 = zeros(1,m);
for k = 1:m
    h_run1(k) = plot(0:T-1, meanA1(:,k), '-',  'LineWidth', 2, 'Color', colors(k,:));
    h_run2(k) = plot(0:T-1, meanA2(:,k), '--', 'LineWidth', 2, 'Color', colors(k,:));
end

% Legenda su 2 colonne:  colonna 1 = tutte le run 1,  colonna 2 = tutte le run 2.
% MATLAB riempie la legenda per RIGHE, quindi si passano gli handle alternati
% run1,run2 per ogni tecnologia (adop k -> stessa riga, run1 a sx, run2 a dx).
leg1 = arrayfun(@(k) ['Adopters [', num2str(k), '], run 1'], 1:m, 'UniformOutput', false);
leg2 = arrayfun(@(k) ['Adopters [', num2str(k), '], run 2'], 1:m, 'UniformOutput', false);
handles = reshape([h_run1; h_run2], 1, []);
labels  = reshape([leg1;  leg2], 1, []);
legend(handles, labels, 'Location', 'best', 'NumColumns', 2, 'Interpreter','latex');
xlabel('Time', 'Interpreter','latex'); ylabel('Fraction of population','Interpreter','latex');
grid on; set(gca,'FontSize',16);
%saveas(gcf, 'caseB_nonunique.png')
