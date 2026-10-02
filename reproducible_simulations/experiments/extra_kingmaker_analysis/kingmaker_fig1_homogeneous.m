% Extra analysis (not in the paper): homogeneous case (n = 1), closed-form
% equilibrium ratio a^[2]*/a^[1]* = (delta^[1]/delta^[2]) (z2/z1) (z1+z3)/(z2+z3),
% with z_k = gamma^[k] x^[k]*, plotted vs z3. Ratio > 1 means technology 3
% reverses the ranking of 1 and 2. Exact for n = 1 (checked against simulate_kingmaker).

clear; close all;

%% Parameters of the incumbents (technologies 1 and 2)
d1 = 0.20;  d2 = 0.30;    % delta^[1], delta^[2]
z1 = 0.20;  z2 = 0.50;    % gamma^[k] x^[k]*

ratio = @(z3) (d1/d2) .* (z2/z1) .* (z1 + z3) ./ (z2 + z3);

%% Limits and rank-reversal threshold
R0    = ratio(0);          % = d1/d2
Rinf  = (d1/d2)*(z2/z1);   % z3 -> inf
z3_star = (z2 - z1*Rinf) / (Rinf - 1);   % ratio(z3_star) = 1

z3 = linspace(0, 1.5, 400);
R  = ratio(z3);

%% Plot
figure; hold on;
yl = [min(0.9,R0*0.95), max(R)*1.05];
if z3_star > 0 && z3_star < max(z3)
    plot([z3_star z3_star], yl, '-', 'LineWidth', 2.0, 'Color', [0.20 0.60 0.25]);
end
plot(z3, R, 'LineWidth', 2.4, 'Color', [0.20 0.30 0.65]);
plot(z3, ones(size(z3)), '--k', 'LineWidth', 1.6);
xlabel('$\gamma^{[3]}\,x^{[3]\star}$', 'Interpreter','latex', 'FontSize',24);
ylabel('$a^{[2]\star}/a^{[1]\star}$', 'Interpreter','latex', 'FontSize',24);
ylim(yl); grid on; set(gca,'FontSize',18);

fprintf('R(0)=%.4f (=delta1/delta2), R(inf)=%.4f, z3* = %.4f\n', R0, Rinf, z3_star);
saveas(gcf, 'kingmaker_fig1_homogeneous.png')
