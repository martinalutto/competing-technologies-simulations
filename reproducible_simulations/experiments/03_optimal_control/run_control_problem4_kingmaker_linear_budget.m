%% ============================================================
%  Problem 4: Kingmaker control after entry of technology 3
%  TUNED FAST VERSION
% ============================================================
% Main speed-ups compared with the full script:
%   1) The pre-entry phase (t=0,...,T1) is simulated only once.
%   2) The forward-backward sweep is run only on the post-entry horizon
%      H = Tfin - T1, where the control can actually affect the system.
%   3) The adjoint recursion uses a direct vector-Jacobian product,
%      avoiding explicit construction of the full Jacobian at each time.
%
% Scenario:
%   - Technologies 1 and 2 are present from t = 0 to T1.
%   - Technology 3 enters at t = T1 with A^{[3]}(T1)=0 and D^{[3]}(T1)=0.
%   - From t = T1 onward, optimize only u^{[3]}.
%   - Objective at Tfin:
%         max  sum_i A_i^{[2]}(Tfin) - sum_i A_i^{[1]}(Tfin)
%     subject to:
%         sum_{t=T1}^{Tfin-2} sum_i u_i^{[3]}(t) <= B,
%         0 <= u_i^{[3]}(t) <= 1.
%
% Kingmaker objective uses only the terminal adoption gap with a linear control budget.

clear; close all; clc;
rng(4);

%% Dimensions and horizons
% Use FAST_TEST = true while tuning parameters. Set it to false to recover
% the full scenario T1=300, Tfin=1000, n=20.
FAST_TEST = false;

if FAST_TEST
    n = 12;
    T1 = 180;
    Tfin = 550;
else
    n = 20;
    T1 = 200;
    Tfin = 400;
end

m = 3;
H = Tfin - T1;       % post-entry horizon
k_new = 3;

%% Network matrices
W = rand(n,n);       W = W ./ sum(W,2);
tilde_W = rand(n,n); tilde_W = tilde_W ./ sum(tilde_W,2);

%% Model parameters
beta  = zeros(n,m);
delta = zeros(n,m);
gamma = zeros(n,m);

% Parameters chosen to make the post-entry FBS better conditioned.
% We keep delta1 < delta2, so technology 1 is more persistent pre-entry.
% Technology 3 has a stronger switching/attraction coefficient gamma3,
% so that it can act as a kingmaker after entry.
beta(:,1)  = 0.45 + 0.20*rand(n,1);
beta(:,2)  = 0.35 + 0.15*rand(n,1);
beta(:,3)  = 0.2 + 0.2*rand(n,1);

delta(:,1) = 0.15 + 0.05*rand(n,1);
delta(:,2) = 0.18 + 0.05*rand(n,1);
delta(:,3) = 0.15 + 0.05*rand(n,1);

gamma(:,1) = 0.25 + 0.04*rand(n,1);
gamma(:,2) = 0.28 + 0.05*rand(n,1);
gamma(:,3) = 0.30 + 0.10*rand(n,1);

% Smaller direct-control coefficient c = 1-lambda-xi avoids very impulsive
% sensitivity of x to u and helps the FBS iterations.
lambda = 0.45*ones(n,m);
xi     = 0.20*ones(n,m);
c      = 1 - lambda - xi;

%% Fixed baseline opinions used in x_U
X0 = zeros(n,m);
X0(:,1) = 0.45 + 0.10*rand(n,1);
X0(:,2) = 0.3 + 0.10*rand(n,1);
X0(:,3) = 0.40 + 0.10*rand(n,1);

%% Initial conditions at t = 0, technologies 1 and 2 only
A0 = zeros(n,m);
D0 = zeros(n,m);
Xinit = X0;

A0(:,1) = 0.090 + 0.015*rand(n,1);
A0(:,2) = 0.1 + 0.010*rand(n,1);
A0(:,3) = 0;
S0 = 1 - sum(A0,2);
S0 = clamp01(S0);
y0 = packState(S0,A0,D0,Xinit,n,m);

%% Control parameters for Problem 4
% Linear budget: sum u_i^{[3]}(t) <= B. With a terminal objective and a
% linear budget, the PMP update is threshold/bang-bang. A moderate budget
% avoids too many active controls and helps convergence.
B = 0.2*n*(H-1);     % about 8% of full post-entry control capacity
umin = 0;
umax = 1;

%% Forward-backward sweep parameters
maxIter = 150;         % tuned for fast convergence during parameter search
tol = 1e-3;
omega = 0.30;        % conservative relaxation, reduces oscillations
useClamping = false; % set true only as a numerical safeguard

%% Parameters structure
p.n = n; p.m = m; p.T1 = T1; p.Tfin = Tfin; p.H = H; p.k_new = k_new;
p.W = W; p.tilde_W = tilde_W;
p.beta = beta; p.delta = delta; p.gamma = gamma;
p.lambda = lambda; p.xi = xi; p.c = c; p.X0 = X0;
p.B = B; p.umin = umin; p.umax = umax;
p.maxIter = maxIter; p.tol = tol; p.omega = omega;
p.useClamping = useClamping;
p.dimY = n + 3*n*m;

%% Pre-entry phase: simulate only once, technologies 1 and 2 active
Y_pre = forwardSimulationPreEntry(y0,p);
yEntry = Y_pre(:,end);

%% Baseline post-entry simulation with no control
Uzero_post = zeros(H,n,m);
Y_un_post = forwardSimulationPost(yEntry,Uzero_post,p);
Y_un = [Y_pre, Y_un_post(:,2:end)];
[S_un,A_un,D_un,X_un] = unpackTrajectory(Y_un,n,m,Tfin);

%% Problem 4: kingmaker control on technology 3 after entry
[Y_ctrl_post,U_post,MU_post,err,objVals,budgetVals] = forwardBackwardSweepKingmakerFast(yEntry,p);
Y_ctrl = [Y_pre, Y_ctrl_post(:,2:end)];
U_ctrl = zeros(Tfin,n,m);
U_ctrl(T1+1:Tfin,:,:) = U_post;
[S_ctrl,A_ctrl,D_ctrl,X_ctrl] = unpackTrajectory(Y_ctrl,n,m,Tfin);

%% Diagnostics
fprintf('\n=== Problem 4: Kingmaker control, fast version ===\n');
fprintf('Pre-entry gap at T1, uncontrolled: mean A2 - mean A1 = %.4f\n', ...
    mean(A_un(T1+1,:,2)) - mean(A_un(T1+1,:,1)));
fprintf('Final gap uncontrolled: mean A2 - mean A1 = %.4f\n', ...
    mean(A_un(Tfin+1,:,2)) - mean(A_un(Tfin+1,:,1)));
fprintf('Final gap controlled:   mean A2 - mean A1 = %.4f\n', ...
    mean(A_ctrl(Tfin+1,:,2)) - mean(A_ctrl(Tfin+1,:,1)));
fprintf('Budget used = %.4f / %.4f\n', computeLinearBudgetPost(U_post,p), B);
fprintf('Max control outside technology 3 = %.4e\n', max(max(max(abs(U_ctrl(:,:,[1,2]))))));
fprintf('Max control before entry = %.4e\n', max(max(max(abs(U_ctrl(1:T1,:,:))))));
fprintf('Max control at final ineffective step u(Tfin-1) = %.4e\n', max(max(abs(U_ctrl(Tfin,:,:)))));
printStateDiagnostics('uncontrolled',S_un,A_un,D_un,X_un);
printStateDiagnostics('controlled',S_ctrl,A_ctrl,D_ctrl,X_ctrl);
printPMPDiagnostics(U_post,p);

%% Plots
plotResultsProblem4(S_un,A_un,D_un,X_un,S_ctrl,A_ctrl,D_ctrl,X_ctrl,U_ctrl,err,budgetVals,p);

%% ============================================================
%  Local functions
% ============================================================

function [Y,U,MU,err,objVals,budgetVals] = forwardBackwardSweepKingmakerFast(yEntry,p)
    H = p.H; n = p.n; m = p.m;
    U = zeros(H,n,m);
    err = zeros(p.maxIter,1);
    objVals = zeros(p.maxIter,1);
    budgetVals = zeros(p.maxIter,1);

    for h = 1:p.maxIter
        Y = forwardSimulationPost(yEntry,U,p);
        MU = backwardAdjointPost(Y,U,p);
        U_tilde = updateControlKingmakerPost(MU,p);

        U_new = p.omega*U_tilde + (1-p.omega)*U;
        U_new = enforceKingmakerAdmissibilityPost(U_new,p);

        err(h) = max(abs(U_new(:)-U(:)));
        objVals(h) = computeObjectiveKingmakerPost(Y,U,p);
        budgetVals(h) = computeLinearBudgetPost(U_new,p);

        if mod(h,10)==0 || h==1
            fprintf('Iter %3d: err = %.3e, J = %.4f, budget = %.2f / %.2f\n', ...
                h,err(h),objVals(h),budgetVals(h),p.B);
        end

        U = U_new;
        if err(h) < p.tol
            err = err(1:h); objVals = objVals(1:h); budgetVals = budgetVals(1:h);
            break;
        end
    end

    Y = forwardSimulationPost(yEntry,U,p);
    MU = backwardAdjointPost(Y,U,p);
end

function Y = forwardSimulationPreEntry(y0,p)
    Y = zeros(p.dimY,p.T1+1);
    Y(:,1) = y0;
    activeTech = 1:2;
    uMat = zeros(p.n,p.m);
    for tIdx = 1:p.T1
        Y(:,tIdx+1) = stateUpdate(Y(:,tIdx),uMat,activeTech,p);
    end
end

function Y = forwardSimulationPost(yEntry,U,p)
    Y = zeros(p.dimY,p.H+1);
    Y(:,1) = yEntry;
    activeTech = 1:3;
    for tIdx = 1:p.H
        uMat = squeeze(U(tIdx,:,:));
        Y(:,tIdx+1) = stateUpdate(Y(:,tIdx),uMat,activeTech,p);
    end
end

function MU = backwardAdjointPost(Y,U,p) %#ok<INUSD>
    % Fast adjoint recursion without explicitly forming the Jacobian.
    % This computes MU(:,t) = J_y(t)' * MU(:,t+1) directly.
    MU = zeros(p.dimY,p.H+1);
    n = p.n;

    % Terminal payoff: sum_i A_i^{[2]}(Tfin) - sum_i A_i^{[1]}(Tfin)
    MU(n + n + (1:n), p.H+1) =  1;   % A^[2]
    MU(n     + (1:n), p.H+1) = -1;   % A^[1]

    for tIdx = p.H:-1:1
        MU(:,tIdx) = adjointStepNoJacobian(Y(:,tIdx),MU(:,tIdx+1),p);
    end
end

function mu = adjointStepNoJacobian(y,muNext,p)
    % Direct vector-Jacobian product for the post-entry dynamics.
    % Direct vector-Jacobian product: avoids building the full Jacobian matrix.
    n = p.n; m = p.m;
    [S,A,D,X] = unpackState(y,n,m);

    muS = muNext(1:n);
    muA = reshape(muNext(n+1:n+n*m),n,m);
    muD = reshape(muNext(n+n*m+1:n+2*n*m),n,m);
    muX = reshape(muNext(n+2*n*m+1:n+3*n*m),n,m);

    gradS = zeros(n,1);
    gradA = zeros(n,m);
    gradD = zeros(n,m);
    gradX = zeros(n,m);

    Wa = p.W*A;
    force = sum(p.beta .* X .* Wa,2);

    % Contributions from S_next = S .* (1 - force)
    gradS = gradS + muS .* (1 - force);
    for k = 1:m
        gradX(:,k) = gradX(:,k) - muS .* S .* p.beta(:,k) .* Wa(:,k);
        gradA(:,k) = gradA(:,k) - p.W' * (muS .* S .* p.beta(:,k) .* X(:,k));
    end

    % Contributions from A_next
    for k = 1:m
        crossD = sum(D,2) - D(:,k);
        tmp = muA(:,k) .* p.beta(:,k) .* X(:,k) .* S;

        gradA(:,k) = gradA(:,k) + muA(:,k).*(1 - p.delta(:,k)) + p.W' * tmp;
        gradS = gradS + muA(:,k) .* p.beta(:,k) .* X(:,k) .* Wa(:,k);
        gradX(:,k) = gradX(:,k) + muA(:,k) .* (p.beta(:,k).*S.*Wa(:,k) + p.gamma(:,k).*crossD);

        for h = 1:m
            if h ~= k
                gradD(:,h) = gradD(:,h) + muA(:,k) .* p.gamma(:,k) .* X(:,k);
            end
        end
    end

    % Contributions from D_next
    totalGX = sum(p.gamma .* X,2);
    for k = 1:m
        switchingOut = totalGX - p.gamma(:,k).*X(:,k);
        gradD(:,k) = gradD(:,k) + muD(:,k).*(1 - switchingOut);
        gradA(:,k) = gradA(:,k) + muD(:,k).*p.delta(:,k);
    end
    for h = 1:m
        tmp = zeros(n,1);
        for k = 1:m
            if k ~= h
                tmp = tmp + muD(:,k).*D(:,k);
            end
        end
        gradX(:,h) = gradX(:,h) - p.gamma(:,h).*tmp;
    end

    % Contributions from X_next
    for k = 1:m
        gradX(:,k) = gradX(:,k) + p.tilde_W' * (muX(:,k).*p.lambda(:,k));
        gradA(:,k) = gradA(:,k) + p.W' * (muX(:,k).*p.xi(:,k));
    end

    mu = packState(gradS,gradA,gradD,gradX,n,m);
end

function U_tilde = updateControlKingmakerPost(MU,p)
    % Linear-budget kingmaker update. Let b_{t,i} be the marginal value of
    % u_i^{[3]}(t). Since the Hamiltonian is affine in u and the budget is
    % linear, the PMP maximizer has a threshold/bang-bang structure: allocate
    % control to the largest positive scores until the budget is exhausted.
    U_tilde = updateControlKingmakerBangBang(MU,p);
    U_tilde = enforceKingmakerAdmissibilityPost(U_tilde,p);
end

function U_tilde = updateControlKingmakerBangBang(MU,p)
    H = p.H; n = p.n; m = p.m; k = p.k_new;
    U_tilde = zeros(H,n,m);

    if p.B <= 0, return; end

    % Score matrix: (H-1) x n.  score(t,i) = c_i^k * (1-X0_i^k) * mu_X_i^k(t+1)
    coeff      = p.c(:,k) .* (1 - p.X0(:,k));           % n x 1
    idxXk      = n + 2*n*m + (k-1)*n + (1:n);           % X^[k] rows in state vector
    scores_mat = (MU(idxXk, 2:H) .* coeff)';            % (H-1) x n

    pos_mask = scores_mat > 0;
    if ~any(pos_mask(:)), return; end

    score_vec           = scores_mat(pos_mask);
    [t_pos, i_pos]      = ind2sub([H-1, n], find(pos_mask));
    [~, ord]            = sort(score_vec, 'descend');

    budgetRemaining = p.B;
    for rr = 1:numel(ord)
        if budgetRemaining <= 0, break; end
        val = min(p.umax, budgetRemaining);
        U_tilde(t_pos(ord(rr)), i_pos(ord(rr)), k) = max(p.umin, val);
        budgetRemaining = budgetRemaining - val;
    end
end

function U = enforceKingmakerAdmissibilityPost(U,p)
    % Post-entry U has rows corresponding to global times T1,...,Tfin-1.
    % Control only technology 3, and not at final ineffective step.
    U(:,:,1) = 0;
    U(:,:,2) = 0;
    U(p.H,:,:) = 0;      % global math t = Tfin-1
    U = min(max(U,p.umin),p.umax);

    Bused = computeLinearBudgetPost(U,p);
    if Bused > p.B && Bused > 0
        U(:,:,p.k_new) = U(:,:,p.k_new) * (p.B/Bused);
    end
end

function J = computeObjectiveKingmakerPost(Y,U,p)
    n = p.n;
    J = sum(Y(n + n + (1:n), p.H+1)) - sum(Y(n + (1:n), p.H+1));
end

function Bused = computeLinearBudgetPost(U,p)
    Bused = sum(U(1:p.H-1, :, p.k_new), 'all');
end

function yNext = stateUpdate(y,u,activeTech,p)
    n = p.n; m = p.m;
    [S,A,D,X] = unpackState(y,n,m);

    % Weighted-average adoption over the network for active technologies
    Wa = zeros(n,m);
    Wa(:,activeTech) = p.W * A(:,activeTech);

    % Susceptible update
    force = sum(p.beta(:,activeTech) .* X(:,activeTech) .* Wa(:,activeTech), 2);
    Snext = S .* (1 - force);

    % A, D, X updates. Inactive technologies keep zeros / X0 (set at init).
    Anext = zeros(n,m);
    Dnext = zeros(n,m);
    Xnext = p.X0;

    totalD  = sum(D(:,activeTech), 2);
    totalGX = sum(p.gamma(:,activeTech) .* X(:,activeTech), 2);

    for k = activeTech
        crossD       = totalD  - D(:,k);
        switchingOut = totalGX - p.gamma(:,k) .* X(:,k);

        Anext(:,k) = A(:,k) + p.beta(:,k).*X(:,k).*S.*Wa(:,k) ...
            - p.delta(:,k).*A(:,k) + p.gamma(:,k).*X(:,k).*crossD;
        Dnext(:,k) = D(:,k) .* (1 - switchingOut) + p.delta(:,k).*A(:,k);

        XU         = p.X0(:,k).*(1 - u(:,k)) + u(:,k);
        Xnext(:,k) = p.c(:,k).*XU + p.lambda(:,k).*(p.tilde_W*X(:,k)) + p.xi(:,k).*Wa(:,k);
    end

    if p.useClamping
        Snext = clamp01(Snext);
        Anext = clamp01(Anext);
        Dnext = clamp01(Dnext);
        Xnext = clamp01(Xnext);
    end

    yNext = packState(Snext,Anext,Dnext,Xnext,n,m);
end


function [S,A,D,X] = unpackState(y,n,m)
    S = y(1:n);
    A = reshape(y(n+1:n+n*m),n,m);
    D = reshape(y(n+n*m+1:n+2*n*m),n,m);
    X = reshape(y(n+2*n*m+1:n+3*n*m),n,m);
end

function y = packState(S,A,D,X,n,m)
    y = [S(:); reshape(A,n*m,1); reshape(D,n*m,1); reshape(X,n*m,1)];
end

function [S,A,D,X] = unpackTrajectory(Y,n,m,~)
    Nt = size(Y,2);   % Tfin+1
    S = Y(1:n, :)';
    A = permute(reshape(Y(n+1       : n+n*m,   :), n, m, Nt), [3,1,2]);
    D = permute(reshape(Y(n+n*m+1   : n+2*n*m, :), n, m, Nt), [3,1,2]);
    X = permute(reshape(Y(n+2*n*m+1 : end,     :), n, m, Nt), [3,1,2]);
end

function plotResultsProblem4(S_un,A_un,D_un,X_un,S_ctrl,A_ctrl,D_ctrl,X_ctrl,U_ctrl,err,budgetVals,p) %#ok<INUSD>
    Tfin = p.Tfin; T1 = p.T1; m = p.m;
    timeState = 0:Tfin;
    timeCtrl = 0:Tfin-1;
    colorsTech = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

    meanA_un = squeeze(mean(A_un,2));
    meanA_ctrl = squeeze(mean(A_ctrl,2));
    meanD_un = squeeze(mean(D_un,2));
    meanD_ctrl = squeeze(mean(D_ctrl,2));
    meanU3 = squeeze(mean(U_ctrl(:,:,3),2));

    figure;
    tiledlayout(2,1);

    nexttile; hold on;
    for k = 1:m
        plot(timeState,meanA_un(:,k),'--','LineWidth',1.5,'Color',colorsTech(k,:));
        plot(timeState,meanA_ctrl(:,k),'-','LineWidth',2.1,'Color',colorsTech(k,:));
    end
    xline(T1,':k','Entry of technology 3','Interpreter','latex','LineWidth',1.3);
    xlabel('Time','Interpreter','latex'); ylabel('Mean adoption','Interpreter','latex');
    title('Problem 4: adoption trajectories','Interpreter','latex');
    grid on; set(gca,'FontSize',14);
    legend({'Unctrl $A^{[1]}$','Ctrl $A^{[1]}$','Unctrl $A^{[2]}$','Ctrl $A^{[2]}$', ...
            'Unctrl $A^{[3]}$','Ctrl $A^{[3]}$'},'Interpreter','latex','Location','bestoutside','NumColumns',2);

    nexttile; hold on;
    gap_un = meanA_un(:,2)-meanA_un(:,1);
    gap_ctrl = meanA_ctrl(:,2)-meanA_ctrl(:,1);
    plot(timeState,gap_un,'--k','LineWidth',1.8);
    plot(timeState,gap_ctrl,'-k','LineWidth',2.4);
    yline(0,':','LineWidth',1.2);
    xline(T1,':k','LineWidth',1.3);
    xlabel('Time','Interpreter','latex'); ylabel('$\overline A^{[2]}-\overline A^{[1]}$','Interpreter','latex');
    title('Aggregate rank gap','Interpreter','latex');
    legend('Uncontrolled','Controlled','Interpreter','latex','Location','best');
    grid on; set(gca,'FontSize',14);

    figure;
    tiledlayout(2,1);

    nexttile; hold on;
    for k = 1:m
        plot(timeState,meanD_un(:,k),'--','LineWidth',1.5,'Color',colorsTech(k,:));
        plot(timeState,meanD_ctrl(:,k),'-','LineWidth',2.1,'Color',colorsTech(k,:));
    end
    xline(T1,':k','LineWidth',1.3);
    xlabel('Time','Interpreter','latex'); ylabel('Mean dissatisfied','Interpreter','latex');
    title('Dissatisfied populations','Interpreter','latex');
    grid on; set(gca,'FontSize',14);

    nexttile; hold on;
    plot(timeCtrl,meanU3,'LineWidth',2.2);
    xline(T1,':k','Entry','Interpreter','latex','LineWidth',1.3);
    xlabel('Time','Interpreter','latex'); ylabel('Mean control on technology 3','Interpreter','latex');
    title('Kingmaker control $u^{[3]}$','Interpreter','latex');
    grid on; set(gca,'FontSize',14);

    figure;
    tiledlayout(2,1);
    nexttile;
    semilogy(err,'LineWidth',2);
    xlabel('Iteration','Interpreter','latex'); ylabel('Control update error','Interpreter','latex');
    title('Forward-backward sweep convergence','Interpreter','latex');
    grid on; set(gca,'FontSize',14);

    nexttile;
    plot(budgetVals,'LineWidth',2); hold on; yline(p.B,'--','LineWidth',1.5);
    xlabel('Iteration','Interpreter','latex'); ylabel('Budget used','Interpreter','latex');
    title('Linear budget usage','Interpreter','latex');
    grid on; set(gca,'FontSize',14);
end

function printStateDiagnostics(name,S,A,D,X)
    fprintf('State bounds (%s): S=[%.3e, %.3e], A=[%.3e, %.3e], D=[%.3e, %.3e], X=[%.3e, %.3e]\n', ...
        name,min(S(:)),max(S(:)),min(A(:)),max(A(:)),min(D(:)),max(D(:)),min(X(:)),max(X(:)));
end

function printPMPDiagnostics(U,p)
    k = p.k_new;
    active = squeeze(U(1:p.H-1,:,k));
    tolBang = 1e-3;
    frac = active(active > tolBang & active < p.umax-tolBang);
    fprintf('PMP check: budget %.4f / %.4f, fractional active controls = %d\n', ...
        computeLinearBudgetPost(U,p),p.B,numel(frac));
end

function val = clamp01(val)
    val = min(max(val,0),1);
end

function idx = indexS(i)
    idx = i;
end

function idx = indexA(i,k,n)
    idx = n + (k-1)*n + i;
end

function idx = indexD(i,k,n,m)
    idx = n + n*m + (k-1)*n + i;
end

function idx = indexX(i,k,n,m)
    idx = n + 2*n*m + (k-1)*n + i;
end
