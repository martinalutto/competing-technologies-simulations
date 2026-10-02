function val = hyp_unique_value(lambda, xi, delta)
% HYP_UNIQUE_VALUE  Left-hand side of the uniqueness condition, eq. (9) of Theorem 2.
%   val = max_k [ max_i xi_i^k / (1 - max_i lambda_i^k)
%                 * max_i (1/delta_i^k) * (1 + 1/(2*min_r delta_i^r)) ]
%   lambda, xi, delta: n x m (agents x technologies). Condition holds iff val < 1.

m = size(lambda, 2);
min_delta_i = min(delta, [], 2);   % min_r delta_i^r

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
