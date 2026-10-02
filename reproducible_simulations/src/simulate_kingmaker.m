function [Aend, Xend, meanA, meanD, meanX, A, D, X, S] = simulate_kingmaker(P)
% SIMULATE_KINGMAKER  Adoption-opinion dynamics with staggered technology entry.
%
% Vectorized over agents. Reproduces the loop in kingmaker.m but with a
% per-technology entry time t_entry(k): technology k is completely absent
% (its a,d,x are frozen at the initial values and it contributes to NO
% cross-technology term) until t >= t_entry(k). This isolates the effect of
% the entering "kingmaker" technology on the post-entry equilibrium.
%
% P is a struct with fields:
%   n, m, T            sizes and horizon
%   t_entry            1 x m entry times (1 = present from the start)
%   beta,delta,gamma   n x m rates
%   lambda,xi          n x m opinion-update weights (lambda+xi < 1)
%   W, tilde_W         n x n row-stochastic matrices
%   A0, X0             n x m initial adopters / initial opinions (x0 anchor)
%   S0                 1 x n initial susceptibles (optional; else 1-sum_k A0)
%   D0                 n x m initial dissatisfied (optional; else 0)
%
% Outputs:
%   Aend  1 x m aggregate adoption A^[k] = sum_i a_i^[k](T)   (as in the notes)
%   Xend  1 x m mean opinion x^[k] at final time
%   meanA,meanD,meanX  T x m node-averaged trajectories
%   A,D,X,S            full state arrays (T x n x m, and T x n for S)

n = P.n; m = P.m; T = P.T;
beta = P.beta; delta = P.delta; gamma = P.gamma;
lambda = P.lambda; xi = P.xi; W = P.W; tilde_W = P.tilde_W;

if isfield(P,'t_entry'); t_entry = P.t_entry; else; t_entry = ones(1,m); end

S = zeros(T,n); A = zeros(T,n,m); D = zeros(T,n,m); X = zeros(T,n,m);
A(1,:,:) = P.A0;
X(1,:,:) = P.X0;
if isfield(P,'D0'); D(1,:,:) = P.D0; else; D(1,:,:) = 0; end
if isfield(P,'S0'); S(1,:) = P.S0; else; S(1,:) = max(0, 1 - sum(A(1,:,:),3)); end

X1 = reshape(X(1,:,:), n, m);   % fixed opinion anchor x0 (n x m, robust for n=1)

for t = 1:T-1
    active = (t >= t_entry);            % 1 x m logical: technologies present at t
    Xt = reshape(X(t,:,:), n, m);
    At = reshape(A(t,:,:), n, m);
    Dt = reshape(D(t,:,:), n, m);
    St = S(t,:)';

    % inactive technologies contribute nothing to any coupling term
    Aeff = At; Aeff(:,~active) = 0;
    Deff = Dt; Deff(:,~active) = 0;
    gX   = gamma .* Xt; gX(:,~active) = 0;

    Wa = W * Aeff;                       % n x m
    Xneighbors = tilde_W * Xt;

    % carry over inactive technologies unchanged
    Xnext = Xt; Anext = At; Dnext = Dt;

    % 1) opinions (active technologies only)
    Xup = (1 - lambda - xi) .* X1 + lambda .* Xneighbors + xi .* Wa;
    Xup = max(0, min(1, Xup));
    Xnext(:,active) = Xup(:,active);

    % 2) susceptibles: force accumulates over active technologies
    adoption_force = sum(beta(:,active) .* Xt(:,active) .* Wa(:,active), 2);
    Snext = St - St .* adoption_force;
    Snext = max(0, min(1, Snext));

    % 3) adopters (active technologies only)
    dissatisfied_others = sum(Deff,2) - Deff;                 % sum_{h~=k} d^[h]
    Aup = At + beta .* Xt .* St .* Wa - delta .* At + gamma .* Xt .* dissatisfied_others;
    Aup = max(0, min(1, Aup));
    Anext(:,active) = Aup(:,active);

    % 4) dissatisfied (active technologies only)
    switching_out = sum(gX,2) - gX;                           % sum_{h~=k} gamma^[h] x^[h]
    Dup = Dt + delta .* At - Dt .* switching_out;
    Dup = max(0, min(1, Dup));
    Dnext(:,active) = Dup(:,active);

    X(t+1,:,:) = Xnext;
    S(t+1,:)   = Snext';
    A(t+1,:,:) = Anext;
    D(t+1,:,:) = Dnext;
end

meanA = squeeze(mean(A,2));
meanD = squeeze(mean(D,2));
meanX = squeeze(mean(X,2));

Aend = squeeze(sum(A(T,:,:),2))';        % aggregate A^[k] = sum_i a_i^[k](T)
Xend = squeeze(mean(X(T,:,:),2))';
end
