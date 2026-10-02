%%%%%%%%%%%%%%%%%%%% KINGMAKER — Figure 2 e 3 (modello eterogeneo su rete)
%% Sweep sui due parametri della tecnologia entrante:
%%   gamma^[3] : attrattivita' di tech 3 per gli utenti dissatisfied
%%   xi^[3]    : intensita' del feedback adoption-to-opinion di tech 3
%%
%% Per ogni coppia (gamma^[3], xi^[3]) si simula la dinamica post-entry fino a
%% convergenza e si calcolano le adozioni aggregate all'equilibrio
%%   A^[k] = sum_i a_i^[k](T).
%%
%% Figura 2 (mappa di fase): heatmap di  DeltaA = A^[2] - A^[1], con contorno
%%   DeltaA = 0 (soglia di rank reversal tra tech 1 e tech 2).
%% Figura 3 (regione kingmaker): mappa binaria
%%   K = 1  se  A^[2] > A^[1]  e  A^[3] < A^[2]   (rank reversal senza che
%%          tech 3 diventi dominante), 0 altrimenti.

clear; close all;
rng(3);

%% ---- Parametri base (tech 1 e tech 2 incumbent, rete eterogenea) ----
n = 30; m = 3; T = 800; t_entry = [1 1 300];

beta = rand(n,m);
beta(:,1) = 0.30 + 0.20*rand(n,1);
beta(:,2) = 0.35 + 0.20*rand(n,1);
beta(:,3) = 0.30 + 0.20*rand(n,1);

delta = rand(n,m);
delta(:,1) = 0.15 + 0.10*rand(n,1);
delta(:,2) = 0.20 + 0.10*rand(n,1);
delta(:,3) = 0.20 + 0.10*rand(n,1);

gamma = zeros(n,m);
gamma(:,1) = 0.30 + 0.10*rand(n,1);
gamma(:,2) = 0.50 + 0.10*rand(n,1);   % gamma^[3] impostato nello sweep

lambda = zeros(n,m);
lambda(:,1) = 0.30; lambda(:,2) = 0.30; lambda(:,3) = 0.30;  % lambda3 fisso -> lambda3+xi3<1

xi = zeros(n,m);
xi(:,1) = 0.20; xi(:,2) = 0.20;       % xi^[3] impostato nello sweep

W = rand(n,n); W = W./sum(W,2);
tilde_W = rand(n,n); tilde_W = tilde_W./sum(tilde_W,2);

A0 = zeros(n,m);
A0(:,1) = 0.30 + 0.20*rand(n,1);
A0(:,2) = 0.30 + 0.20*rand(n,1);
X0 = 0.4*ones(n,m);

P.n=n; P.m=m; P.T=T; P.t_entry=t_entry;
P.beta=beta; P.delta=delta; P.lambda=lambda;
P.W=W; P.tilde_W=tilde_W; P.A0=A0; P.X0=X0;

%% ---- Griglia (gamma^[3], xi^[3]) ----
Ng = 60;  Nx = 60;
g3_vec = linspace(0.0, 1.2, Ng);       % asse x: gamma^[3]
x3_vec = linspace(0.0, 0.6, Nx);       % asse y: xi^[3]  (lambda3+xi3 <= 0.9 < 1)

A1g = zeros(Nx, Ng);   % righe = xi3 (y), colonne = gamma3 (x)
A2g = zeros(Nx, Ng);
A3g = zeros(Nx, Ng);

%% ---- Baseline CONTROFATTUALE: nessuna tecnologia 3 ----
% Serve per dire che il reversal e' INDOTTO da tech 3: si stabilisce chi
% vince tra 1 e 2 quando la terza tecnologia non entra mai.
Pbase = P; Pbase.t_entry = [1 1 T+1];   % tech 3 non entra
Pbase.gamma = gamma; Pbase.gamma(:,3) = 0;
Pbase.xi    = xi;    Pbase.xi(:,3)    = 0;
Abase = simulate_kingmaker(Pbase);
A1_base = Abase(1); A2_base = Abase(2);
fprintf('Baseline senza tech 3: A1=%.3f, A2=%.3f -> vince tech %d\n', ...
    A1_base, A2_base, 1 + (A2_base > A1_base));
assert(A1_base > A2_base, ...
    'Baseline non ha tech 1 in testa: il "reversal indotto" non e'' ben definito.');

fprintf('Sweep %d x %d = %d simulazioni...\n', Nx, Ng, Nx*Ng);
tic;
for ix = 1:Nx
    for ig = 1:Ng
        P.gamma = gamma; P.gamma(:,3) = g3_vec(ig);
        P.xi    = xi;    P.xi(:,3)    = x3_vec(ix);
        Aend = simulate_kingmaker(P);
        A1g(ix,ig) = Aend(1);
        A2g(ix,ig) = Aend(2);
        A3g(ix,ig) = Aend(3);
    end
end
fprintf('fatto in %.1f s\n', toc);

DeltaA = A2g - A1g;                         % Figura 2

% Figura 3: classificazione dell'esito rispetto al baseline (A1_base > A2_base).
%   reversal indotto da tech 3  <=>  A2 > A1  (il ranking 1>2 del baseline si ribalta)
%   tech 3 dominante            <=>  A3 >= A2 (la terza scavalca la vincitrice)
reversal = (A2g > A1g);
dominant = (A3g >= A2g);
C = zeros(Nx, Ng);                 % 0 = nessun reversal (tech 1 resta in testa)
C(reversal & ~dominant) = 1;       % 1 = KINGMAKER: 3 induce il reversal ma non domina
C(reversal &  dominant) = 2;       % 2 = reversal ma 3 diventa dominante
K = double(C == 1);                % indicatore kingmaker (come da definizione)

% colormap divergente con i colori dei problemi di ottimizzazione:
%   tech 1 (blu) quando vince la 1, tech 2 (arancio) quando vince la 2.
nc = 256; tt = linspace(0,1,nc)';
c_lo  = [0.00 0.45 0.74];    % tech 1 - blu
c_mid = [0.97 0.97 0.97];    % quasi bianco
c_hi  = [0.85 0.33 0.10];    % tech 2 - arancio
cmapBR = zeros(nc,3);
lo = tt < 0.5;
f = tt(lo)/0.5;              cmapBR(lo,:)  = (1-f).*c_lo  + f.*c_mid;   % blu -> bianco
f = (tt(~lo)-0.5)/0.5;       cmapBR(~lo,:) = (1-f).*c_mid + f.*c_hi;    % bianco -> rosso

%% ================= FIGURA 2: mappa di fase =================
figure; hold on;
imagesc(g3_vec, x3_vec, DeltaA); axis xy; axis tight;
% colormap divergente centrata su 0
cmax = max(abs(DeltaA(:)));
caxis([-cmax cmax]);
colormap(gca, cmapBR);
cb = colorbar; ylabel(cb, '$\Delta A = A^{[2]}-A^{[1]}$', 'Interpreter','latex', 'FontSize',22);
% contorno DeltaA = 0 (soglia di rank reversal)
contour(g3_vec, x3_vec, DeltaA, [0 0], 'k', 'LineWidth', 2.2);
xlabel('$\gamma^{[3]}$', 'Interpreter','latex', 'FontSize',26);
ylabel('$\xi^{[3]}$', 'Interpreter','latex', 'FontSize',26);
set(gca,'FontSize',18);
%saveas(gcf, 'kingmaker_fig2_phasemap.png')

%% ================= FIGURA 3: regione kingmaker (3 categorie) =================
% Distingue i tre esiti dell'ingresso di tech 3 rispetto al baseline (tech 1
% davanti). La banda kingmaker (blu) e' racchiusa tra "niente reversal" e
% "tech 3 dominante": e' esattamente dove 3 INDUCE il reversal SENZA dominare.
figure; hold on;
imagesc(g3_vec, x3_vec, C); axis xy; axis tight;
% colore = tecnologia in testa (palette dei problemi di ottimizzazione)
cmap3 = [0.62 0.80 0.93;    % 0 = niente reversal  -> tech 1 in testa (blu)
         0.96 0.72 0.55;    % 1 = kingmaker        -> tech 2 in testa (arancio)
         0.72 0.86 0.62];   % 2 = tech 3 dominante -> tech 3 in testa (verde)
colormap(gca, cmap3); caxis([-0.5 2.5]);
% bordi: rank reversal (A2=A1) e dominanza di tech 3 (A3=A2)
contour(g3_vec, x3_vec, DeltaA,   [0 0], 'k', 'LineWidth', 2.0);
contour(g3_vec, x3_vec, A3g-A2g,  [0 0], 'k', 'LineWidth', 2.0);
xlabel('$\gamma^{[3]}$', 'Interpreter','latex', 'FontSize',26);
ylabel('$\xi^{[3]}$', 'Interpreter','latex', 'FontSize',26);

textBlue   = [0.00 0.35 0.70];
textOrange = [0.85 0.25 0.05];
textGreen  = [0.25 0.60 0.05];

text(0.09, 0.30, '$A^{[2]}<A^{[1]}$', 'Color',textBlue,'Interpreter','latex', ...
     'FontSize',17, 'Rotation',0, 'HorizontalAlignment','center');
text(0.39, 0.30, '$A^{[2]}>A^{[1]},\ A^{[3]}<A^{[2]}$', 'Interpreter','latex', ...
     'FontSize',17, 'Rotation',0, 'Color',textOrange, 'HorizontalAlignment','center');
text(0.88, 0.30, '$A^{[2]}>A^{[1]},\ A^{[3]}>A^{[2]}$', 'Interpreter','latex', ...
     'FontSize',17, 'Rotation',0, 'Color',textGreen,'HorizontalAlignment','center');
set(gca,'FontSize',18);
saveas(gcf, 'kingmaker_fig3_region.png')

fprintf('Esiti: no-reversal %.1f%% | kingmaker %.1f%% | 3-dominante %.1f%%\n', ...
    100*nnz(C==0)/numel(C), 100*nnz(C==1)/numel(C), 100*nnz(C==2)/numel(C));
