%%%%%%%%%%%%%%%%%%%% ADOZIONE multi-tecnologia — Caso A
%% Regime in cui l'ipotesi di unicita' (eq:hyp-unique) e' SODDISFATTA

rng(1);
n = 30;    % numero di agenti
m = 5;     % numero di tecnologie
T = 350;   % orizzonte temporale

%% Parametri di rete e budget (Assumption (i), (iv))
[W, tilde_W, ~, ~] = build_params_assumption(n, m);
beta  = rand(n,m);   beta  = beta  ./ max(1, sum(beta,2));   % Assumption (iv): sum_k beta_i^k  <= 1
gamma = rand(n,m);   gamma = gamma ./ max(1, sum(gamma,2));  % Assumption (iv): sum_k gamma_i^k <= 1

%% lambda, xi, delta calibrati per soddisfare eq:hyp-unique GENUINAMENTE.
%% Struttura del LHS:  max_k [ max_i xi/(1-max_i lambda) * max_i (1/delta)(1+1/min_r delta) ].
%% Il fattore term2 = max_i (1/delta_i^k)(1+1/min_r delta_i^r) e' >= 2 e cresce
%% come 1/delta^2: per delta piccolo (~0.3) vale ~20, e servirebbe xi~0 (degenere).
%% Per soddisfare la condizione con xi NON trascurabile bisogna quindi tenere
%% delta VICINO A 1 (delta=1 -> term2=2), cosi' term1 = max_i xi/(1-max_i lambda)
%% puo' restare moderato. NB: delta grande = churn elevato -> equilibrio piu' basso;
%% e' la tensione strutturale intrinseca alla condizione di unicita'.
lambda = 0.30 * rand(n,m);                     % lambda_i^k < 1  (Assumption (ii): esiste sempre j con lambda_j^k<1)
delta  = 0.45 + 0.3 * rand(n,m);              
xi     = 0.1 * rand(n,m) .* (1 - lambda);     % xi moderato, lambda_i^k + xi_i^k < 1 per costruzione

val = hyp_unique_value(lambda, xi, delta);
fprintf('[Caso A] eq:hyp-unique LHS = %.4f  (< 1 richiesto)\n', val);
assert(val < 1, 'eq:hyp-unique NON soddisfatta (LHS = %.4f >= 1)', val);

%% Opinioni iniziali strettamente positive (Assumption (iii)), fissate uguali nei due run
X0 = rand(n,m);

%% Due diverse condizioni iniziali di ADOZIONE (non delle opinioni)
A0_run1 = zeros(n,m);
A0_run1(:,1) = 0.1*rand(n,1);
A0_run1(:,2) = 0.05*rand(n,1);
A0_run1(:,5) = 0.1*rand(n,1);

A0_run2 = zeros(n,m);
A0_run2(:,3) = 0.15*rand(n,1);
A0_run2(:,4) = 0.1*rand(n,1);
A0_run2(:,1) = 0.05*rand(n,1);
A0_run2(:,5) = 0.3*rand(n,1);

%% Simulazioni
[~, A1, D1, X1] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0_run1, X0);
[~, A2, ~,  X2] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0_run2, X0);

meanA1 = squeeze(mean(A1,2));
meanA2 = squeeze(mean(A2,2));
meanD1 = squeeze(mean(D1,2));
meanX1 = squeeze(mean(X1,2));
meanX2 = squeeze(mean(X2,2));

fprintf('[Caso A] Frazioni finali adozione (run 1): '); disp(meanA1(end,:));
fprintf('[Caso A] Frazioni finali adozione (run 2): '); disp(meanA2(end,:));
fprintf('[Caso A] Differenza massima equilibrio adozione: %.2e\n', max(abs(meanA1(end,:) - meanA2(end,:))));
fprintf('[Caso A] Differenza massima equilibrio opinioni : %.2e\n', max(abs(meanX1(end,:) - meanX2(end,:))));

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
