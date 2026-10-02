%%%%%%%%%%%%%%%%%%%% KINGMAKER — Figura 1: caso omogeneo (n = 1)
%% Rapporto di adozione all'equilibrio tra tech 2 e tech 1, in forma chiusa:
%%
%%   a^[2]* / a^[1]* = (delta^[1]/delta^[2]) * (z2/z1) * (z1+z3)/(z2+z3),
%%
%% con z_k = gamma^[k] x^[k]* (attrattivita' effettiva della tecnologia k nei
%% flussi di switching all'equilibrio). La figura mostra il rapporto al variare
%% di z3 = gamma^[3] x^[3]*, con la soglia orizzontale a 1: quando il rapporto
%% supera 1 si ha il rank reversal tra le tecnologie 1 e 2 indotto dalla terza.
%%
%% NB: la formula e' esatta nel caso omogeneo (verificata contro la simulazione
%% a singolo nodo di simulate_kingmaker.m a ~1e-16, con S* -> 0).

clear; close all;

%% Parametri di riferimento (tech 1 e tech 2 incumbent)
d1 = 0.20;  d2 = 0.30;    % delta^[1], delta^[2]
z1 = 0.20;  z2 = 0.50;    % z1 = gamma^[1] x^[1]*, z2 = gamma^[2] x^[2]*

ratio = @(z3) (d1/d2) .* (z2/z1) .* (z1 + z3) ./ (z2 + z3);

%% Limiti notevoli
R0    = ratio(0);          % z3 -> 0   : nessuna terza tecnologia  -> = d1/d2
Rinf  = (d1/d2)*(z2/z1);   % z3 -> inf : = (d1/d2)(z2/z1)
% Soglia di rank reversal: ratio(z3*) = 1
z3_star = (z2 - z1*Rinf) / (Rinf - 1);   % risolve (d1/d2)(z2/z1)(z1+z3)/(z2+z3)=1

z3 = linspace(0, 1.5, 400);
R  = ratio(z3);

%% Plot
figure; hold on;
yl = [min(0.9,R0*0.95), max(R)*1.05];
% riga verticale verde in corrispondenza della soglia di rank reversal
if z3_star > 0 && z3_star < max(z3)
    plot([z3_star z3_star], yl, '-', 'LineWidth', 2.0, 'Color', [0.20 0.60 0.25]);
end
plot(z3, R, 'LineWidth', 2.4, 'Color', [0.20 0.30 0.65]);
plot(z3, ones(size(z3)), '--k', 'LineWidth', 1.6);           % soglia a 1
xlabel('$\gamma^{[3]}\,x^{[3]\star}$', 'Interpreter','latex', 'FontSize',24);
ylabel('$a^{[2]\star}/a^{[1]\star}$', 'Interpreter','latex', 'FontSize',24);
ylim(yl); grid on; set(gca,'FontSize',18);

fprintf('R(0)=%.4f (=delta1/delta2), R(inf)=%.4f, z3* = %.4f\n', R0, Rinf, z3_star);
saveas(gcf, 'kingmaker_fig1_homogeneous.png')
