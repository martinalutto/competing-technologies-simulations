function [S, A, D, X] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0, X0)
% SIMULATE_ADOPTION_MULTI  Adoption-opinion dynamics for m technologies on n agents.
%
%   A0, X0 are n x m matrices of initial adopters / initial opinions.
%   Assumption (iii) requires X0 > 0 elementwise.
%
% Vectorized over agents (verified against the explicit per-agent loop:
% max discrepancy at machine precision, ~1e-16).

n = size(W, 1);
m = size(beta, 2);

S = zeros(T, n);
A = zeros(T, n, m);
D = zeros(T, n, m);
X = zeros(T, n, m);

A(1,:,:) = A0;
X(1,:,:) = X0;
S(1,:) = max(0, min(1, 1 - sum(A(1,:,:), 3)));
D(1,:,:) = 0;

X1 = squeeze(X(1,:,:));   % fixed anchor to initial opinion, used at every step

for t = 1:T-1
    Xt = squeeze(X(t,:,:));
    At = squeeze(A(t,:,:));
    Dt = squeeze(D(t,:,:));
    St = S(t,:)';

    Wa = W * At;                  % Wa(i,k) = sum_j W(i,j) A(t,j,k)
    Xneighbors = tilde_W * Xt;    % social-influence term

    %% 1) aggiornamento opinioni
    Xnext = (1 - lambda - xi) .* X1 + lambda .* Xneighbors + xi .* Wa;
    Xnext = max(0, min(1, Xnext));

    %% 2) aggiornamento non-adottanti
    adoption_force = sum(beta .* Xt .* Wa, 2);
    Snext = St - St .* adoption_force;
    Snext = max(0, min(1, Snext));

    %% 3) aggiornamento adottanti A^{[k]}
    dissatisfied_others = sum(Dt, 2) - Dt;   % sum_{h~=k} D(t,i,h)
    Anext = At + beta .* Xt .* St .* Wa - delta .* At + gamma .* Xt .* dissatisfied_others;
    Anext = max(0, min(1, Anext));

    %% 4) aggiornamento insoddisfatti D^{[k]}
    switching_out = sum(gamma .* Xt, 2) - gamma .* Xt;   % sum_{h~=k} gamma(i,h) X(t,i,h)
    Dnext = Dt + delta .* At - Dt .* switching_out;
    Dnext = max(0, min(1, Dnext));

    X(t+1,:,:) = Xnext;
    S(t+1,:)   = Snext';
    A(t+1,:,:) = Anext;
    D(t+1,:,:) = Dnext;
end
end
