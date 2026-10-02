function [S, A, D, X] = simulate_adoption_multi(T, W, tilde_W, beta, gamma, delta, lambda, xi, A0, X0)
% SIMULATE_ADOPTION_MULTI  MATE dynamics, eq. (1a)-(1d), for m technologies on n agents.
%   A0, X0: n x m initial adopters and predispositions p (Assumption 1(ii): X0 > 0).
%   S: T x n susceptibles; A, D, X: T x n x m adopters, dissatisfied, opinions.

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

X1 = squeeze(X(1,:,:));   % predisposition p (opinion anchor)

for t = 1:T-1
    Xt = squeeze(X(t,:,:));
    At = squeeze(A(t,:,:));
    Dt = squeeze(D(t,:,:));
    St = S(t,:)';

    Wa = W * At;
    Xneighbors = tilde_W * Xt;

    % Opinions
    Xnext = (1 - lambda - xi) .* X1 + lambda .* Xneighbors + xi .* Wa;
    Xnext = max(0, min(1, Xnext));

    % Susceptibles
    adoption_force = sum(beta .* Xt .* Wa, 2);
    Snext = St - St .* adoption_force;
    Snext = max(0, min(1, Snext));

    % Adopters
    dissatisfied_others = sum(Dt, 2) - Dt;   % sum_{h~=k} D(t,i,h)
    Anext = At + beta .* Xt .* St .* Wa - delta .* At + gamma .* Xt .* dissatisfied_others;
    Anext = max(0, min(1, Anext));

    % Dissatisfied
    switching_out = sum(gamma .* Xt, 2) - gamma .* Xt;   % sum_{h~=k} gamma(i,h) X(t,i,h)
    Dnext = Dt + delta .* At - Dt .* switching_out;
    Dnext = max(0, min(1, Dnext));

    X(t+1,:,:) = Xnext;
    S(t+1,:)   = Snext';
    A(t+1,:,:) = Anext;
    D(t+1,:,:) = Dnext;
end
end
