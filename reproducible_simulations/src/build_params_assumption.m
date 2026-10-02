function [W, tilde_W, beta, gamma] = build_params_assumption(n, m)
% BUILD_PARAMS_ASSUMPTION  Network and budget parameters satisfying
% Assumption (i) and (iv).
%
%   (i)  W, tilde_W row-stochastic; W irreducible (strongly connected).
%        rand(n,n) has strictly positive entries, so after row-normalizing
%        the associated digraph is complete -> trivially strongly connected.
%  (iv)  sum_k beta_i^k <= 1 and sum_k gamma_i^k <= 1 for every agent i.

W = rand(n, n);
W = W ./ sum(W, 2);

tilde_W = rand(n, n);
tilde_W = tilde_W ./ sum(tilde_W, 2);

beta = rand(n, m);
beta = beta ./ max(1, sum(beta, 2));     % only rescales rows whose sum > 1

gamma = 0.5 * rand(n, m);
gamma = gamma ./ max(1, sum(gamma, 2));  % already <=1 in expectation, kept for safety
end
