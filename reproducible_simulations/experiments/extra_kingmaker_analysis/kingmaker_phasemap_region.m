% Extra analysis (not in the paper): heterogeneous network, technology 3 enters
% at t = 300. Sweep over (gamma^[3], xi^[3]) of the aggregate equilibrium adoption
% A^[k] = sum_i a_i^[k](T). Plots: (1) phase map of A^[2]-A^[1]; (2) kingmaker
% region (A^[2] > A^[1] and A^[3] < A^[2]).

clear; close all;
rng(3);

%% Base parameters (technologies 1 and 2 are incumbents)
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
gamma(:,2) = 0.50 + 0.10*rand(n,1);   % gamma^[3] set in the sweep

lambda = zeros(n,m);
lambda(:,1) = 0.30; lambda(:,2) = 0.30; lambda(:,3) = 0.30;

xi = zeros(n,m);
xi(:,1) = 0.20; xi(:,2) = 0.20;       % xi^[3] set in the sweep

W = rand(n,n); W = W./sum(W,2);
tilde_W = rand(n,n); tilde_W = tilde_W./sum(tilde_W,2);

A0 = zeros(n,m);
A0(:,1) = 0.30 + 0.20*rand(n,1);
A0(:,2) = 0.30 + 0.20*rand(n,1);
X0 = 0.4*ones(n,m);

P.n=n; P.m=m; P.T=T; P.t_entry=t_entry;
P.beta=beta; P.delta=delta; P.lambda=lambda;
P.W=W; P.tilde_W=tilde_W; P.A0=A0; P.X0=X0;

%% Grid (gamma^[3], xi^[3])
Ng = 60;  Nx = 60;
g3_vec = linspace(0.0, 1.2, Ng);
x3_vec = linspace(0.0, 0.6, Nx);       % lambda3 + xi3 < 1

A1g = zeros(Nx, Ng);   % rows = xi3, columns = gamma3
A2g = zeros(Nx, Ng);
A3g = zeros(Nx, Ng);

%% Counterfactual baseline without technology 3 (technology 1 must lead)
Pbase = P; Pbase.t_entry = [1 1 T+1];   % technology 3 never enters
Pbase.gamma = gamma; Pbase.gamma(:,3) = 0;
Pbase.xi    = xi;    Pbase.xi(:,3)    = 0;
Abase = simulate_kingmaker(Pbase);
A1_base = Abase(1); A2_base = Abase(2);
fprintf('Baseline without tech 3: A1=%.3f, A2=%.3f -> tech %d wins\n', ...
    A1_base, A2_base, 1 + (A2_base > A1_base));
assert(A1_base > A2_base, ...
    'Tech 1 does not lead in the baseline: "induced reversal" is not well defined.');

fprintf('Sweep %d x %d = %d simulations...\n', Nx, Ng, Nx*Ng);
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
fprintf('done in %.1f s\n', toc);

DeltaA = A2g - A1g;

% Outcome classes w.r.t. the baseline: 0 = no reversal, 1 = kingmaker, 2 = technology 3 dominant
reversal = (A2g > A1g);
dominant = (A3g >= A2g);
C = zeros(Nx, Ng);
C(reversal & ~dominant) = 1;
C(reversal &  dominant) = 2;
K = double(C == 1);                % kingmaker indicator

% Diverging colormap: blue (technology 1 leads) -> white -> orange (technology 2 leads)
nc = 256; tt = linspace(0,1,nc)';
c_lo  = [0.00 0.45 0.74];
c_mid = [0.97 0.97 0.97];
c_hi  = [0.85 0.33 0.10];
cmapBR = zeros(nc,3);
lo = tt < 0.5;
f = tt(lo)/0.5;              cmapBR(lo,:)  = (1-f).*c_lo  + f.*c_mid;
f = (tt(~lo)-0.5)/0.5;       cmapBR(~lo,:) = (1-f).*c_mid + f.*c_hi;

%% Phase map of DeltaA, with rank-reversal contour DeltaA = 0
figure; hold on;
imagesc(g3_vec, x3_vec, DeltaA); axis xy; axis tight;
cmax = max(abs(DeltaA(:)));
caxis([-cmax cmax]);
colormap(gca, cmapBR);
cb = colorbar; ylabel(cb, '$\Delta A = A^{[2]}-A^{[1]}$', 'Interpreter','latex', 'FontSize',22);
contour(g3_vec, x3_vec, DeltaA, [0 0], 'k', 'LineWidth', 2.2);
xlabel('$\gamma^{[3]}$', 'Interpreter','latex', 'FontSize',26);
ylabel('$\xi^{[3]}$', 'Interpreter','latex', 'FontSize',26);
set(gca,'FontSize',18);

%% Kingmaker region (3 outcome classes)
figure; hold on;
imagesc(g3_vec, x3_vec, C); axis xy; axis tight;
cmap3 = [0.62 0.80 0.93;    % 0: no reversal (technology 1 leads)
         0.96 0.72 0.55;    % 1: kingmaker (technology 2 leads)
         0.72 0.86 0.62];   % 2: technology 3 dominant
colormap(gca, cmap3); caxis([-0.5 2.5]);
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

fprintf('Outcomes: no reversal %.1f%% | kingmaker %.1f%% | 3 dominant %.1f%%\n', ...
    100*nnz(C==0)/numel(C), 100*nnz(C==1)/numel(C), 100*nnz(C==2)/numel(C));
