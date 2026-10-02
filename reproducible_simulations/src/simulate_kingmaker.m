function [Aend, Xend, meanA, meanD, meanX, A, D, X, S] = simulate_kingmaker(P)
% SIMULATE_KINGMAKER  MATE dynamics with staggered technology entry.
%   Technology k is frozen at its initial state and enters no coupling term
%   until t >= P.t_entry(k).
%
% P fields: n, m, T; t_entry (1 x m, optional, default all 1); beta, delta,
%   gamma, lambda, xi (n x m); W, tilde_W (n x n); A0, X0 (n x m, X0 = p);
%   S0 (1 x n, optional); D0 (n x m, optional).
% Outputs: Aend, Xend (1 x m) aggregate adoption sum_i a_i^[k](T) and mean
%   final opinion; meanA, meanD, meanX (T x m) node averages; A, D, X, S full states.

n = P.n; m = P.m; T = P.T;
beta = P.beta; delta = P.delta; gamma = P.gamma;
lambda = P.lambda; xi = P.xi; W = P.W; tilde_W = P.tilde_W;

if isfield(P,'t_entry'); t_entry = P.t_entry; else; t_entry = ones(1,m); end

S = zeros(T,n); A = zeros(T,n,m); D = zeros(T,n,m); X = zeros(T,n,m);
A(1,:,:) = P.A0;
X(1,:,:) = P.X0;
if isfield(P,'D0'); D(1,:,:) = P.D0; else; D(1,:,:) = 0; end
if isfield(P,'S0'); S(1,:) = P.S0; else; S(1,:) = max(0, 1 - sum(A(1,:,:),3)); end

X1 = reshape(X(1,:,:), n, m);   % predisposition p (reshape is safe for n = 1)

for t = 1:T-1
    active = (t >= t_entry);
    Xt = reshape(X(t,:,:), n, m);
    At = reshape(A(t,:,:), n, m);
    Dt = reshape(D(t,:,:), n, m);
    St = S(t,:)';

    % inactive technologies do not enter the coupling terms
    Aeff = At; Aeff(:,~active) = 0;
    Deff = Dt; Deff(:,~active) = 0;
    gX   = gamma .* Xt; gX(:,~active) = 0;

    Wa = W * Aeff;
    Xneighbors = tilde_W * Xt;

    Xnext = Xt; Anext = At; Dnext = Dt;

    % Opinions
    Xup = (1 - lambda - xi) .* X1 + lambda .* Xneighbors + xi .* Wa;
    Xup = max(0, min(1, Xup));
    Xnext(:,active) = Xup(:,active);

    % Susceptibles
    adoption_force = sum(beta(:,active) .* Xt(:,active) .* Wa(:,active), 2);
    Snext = St - St .* adoption_force;
    Snext = max(0, min(1, Snext));

    % Adopters
    dissatisfied_others = sum(Deff,2) - Deff;
    Aup = At + beta .* Xt .* St .* Wa - delta .* At + gamma .* Xt .* dissatisfied_others;
    Aup = max(0, min(1, Aup));
    Anext(:,active) = Aup(:,active);

    % Dissatisfied
    switching_out = sum(gX,2) - gX;
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

Aend = squeeze(sum(A(T,:,:),2))';
Xend = squeeze(mean(X(T,:,:),2))';
end
