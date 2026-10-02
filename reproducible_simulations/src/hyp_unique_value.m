function val = hyp_unique_value(lambda, xi, delta)
% HYP_UNIQUE_VALUE  LHS of the uniqueness condition (eq:hyp-unique).
%
%   val = max_k ( max_i xi_i^k / (1 - max_i lambda_i^k) ...
%                 * max_i [ (1/delta_i^k) * (1 + 1/min_r delta_i^r) ] )
%
% lambda, xi, delta are n x m matrices (rows = agents, cols = technologies).
% val < 1  =>  uniqueness hypothesis satisfied.

m = size(lambda, 2);
min_delta_i = min(delta, [], 2);   % n x 1, min_r delta_i^r for each agent i

vals_k = zeros(1, m);
for k = 1:m
    max_lambda_k = max(lambda(:,k));
    max_xi_k     = max(xi(:,k));
    term1 = max_xi_k / (1 - max_lambda_k);
    term2 = max( (1 ./ delta(:,k)) .* (1 + 1 ./ (2*min_delta_i)) );
    vals_k(k) = term1 * term2;
end
val = max(vals_k);
end
