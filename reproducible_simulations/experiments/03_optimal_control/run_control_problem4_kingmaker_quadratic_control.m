%% ============================================================
%  Problem 4: Kingmaker control after entry of technology 3
%  QUADRATIC-COST CONTROL VERSION
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
%     with quadratic intervention cost:
%         max gap(Tfin) - qK * sum_{t=T1}^{Tfin-2} sum_i (u_i^{[3]}(t))^2,
%         0 <= u_i^{[3]}(t) <= 1.
%
% Kingmaker objective uses the terminal adoption gap with a quadratic intervention cost.

clear; close all; clc;
rng(2);

%% Dimensions and horizons
% Use FAST_TEST = true while tuning parameters.
FAST_TEST = false;

if FAST_TEST
    n = 12;
    T1 = 30;
    Tfin = 250;
else
    n = 20;
    T1 = 60;    % short pre-entry: tech 2 still leads at T1 (gap > 0) before gamma^[1] erodes it
    Tfin = 150; % long post-entry (H=380): gap flips negative without control, positive with it
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

% Kingmaker scenario (true flip):
%   WITHOUT control: the large D^[2] pool slowly refills A^[1] post-T1
%                    (via gamma^[1]*X^[1]*D^[2]) and the gap A2-A1 flips negative.
%   WITH control:    tech 3 acts as a "sponge" for D^[2], absorbing it before
%                    it reaches A^[1]; D^[3] then refills A^[2] preferentially
%                    (gamma^[2]*X^[2] >> gamma^[1]*X^[1]) and the gap stays positive.
%   X0^[3] is kept near zero so that uncontrolled tech 3 stays dormant (X^[3]≈0.07);
%   only the control u pushes X^[3] toward 1, activating the sponge mechanism.
beta(:,1)  = 0.20 + 0.08*rand(n,1);   % challenger: low adoption spread
beta(:,2)  = 0.55 + 0.10*rand(n,1);   % incumbent: high adoption spread
beta(:,3)  = 0.25 + 0.10*rand(n,1);   % new entrant

delta(:,1) = 0.08 + 0.03*rand(n,1);   % challenger: sticky (low churn)
delta(:,2) = 0.22 + 0.05*rand(n,1);   % incumbent: leaky (high churn → large D^[2])
delta(:,3) = 0.15 + 0.05*rand(n,1);

% gamma^[1] moderate: large enough to flip the gap post-T1 without control,
% small enough that short T1=50 steps haven't let tech 1 overtake tech 2 yet.
gamma(:,1) = 0.1 + 0.04*rand(n,1);   % challenger: moderate attraction to switchers
gamma(:,2) = 0.30 + 0.05*rand(n,1);   % incumbent: good at re-attracting
gamma(:,3) = 0.5 + 0.10*rand(n,1);   % kingmaker: strong attractor

% Smaller direct-control coefficient c = 1-lambda-xi avoids very impulsive
% sensitivity of x to u and helps the FBS iterations.
lambda = 0.3*ones(n,m);
xi     = 0.20*ones(n,m);
c      = 1 - lambda - xi;

%% Fixed baseline opinions used in x_U
X0 = zeros(n,m);
X0(:,1) = 0.28 + 0.08*rand(n,1);   % challenger
X0(:,2) = 0.50 + 0.08*rand(n,1);   % incumbent
X0(:,3) = 0.05 + 0.05*rand(n,1);   % near-zero: tech 3 dormant without control
                                     % (X^[3]_uncontrolled ≈ 0.07; with u=1 → X^[3] ≈ 0.65)

%% Initial conditions at t = 0, technologies 1 and 2 only
A0 = zeros(n,m);
D0 = zeros(n,m);
Xinit = X0;

A0(:,1) = 0.02  + 0.01*rand(n,1);   % tech 1: challenger (tiny initial share)
A0(:,2) = 0.2  + 0.03*rand(n,1);   % tech 2: incumbent (large initial share)
A0(:,3) = 0;
S0 = 1 - sum(A0,2);
S0 = clamp01(S0);
y0 = packState(S0,A0,D0,Xinit,n,m);

%% Control parameters for Problem 4
% Quadratic running cost on the kingmaker intervention:
%   J = final gap - qK * sum_t sum_i (u_i^{[3]}(t))^2.
% This gives a continuous projected PMP update instead of a bang-bang one.
qK = 0.01;           % quadratic penalty: u* = b/(2*qK). Smaller → stronger controls.
                     % 0.1 → u_avg ≈ 1-2% (invisible); 0.01 → u_avg ≈ 10-30% (visible)
umin = 0;
umax = 1;

%% Forward-backward sweep parameters
maxIter = 100;         % tuned for fast convergence during parameter search
tol = 1e-5;
omega = 0.20;        % conservative relaxation, reduces oscillations
useClamping = false; % set true only as a numerical safeguard

%% Parameters structure
p.n = n; p.m = m; p.T1 = T1; p.Tfin = Tfin; p.H = H; p.k_new = k_new;
p.W = W; p.tilde_W = tilde_W;
p.beta = beta; p.delta = delta; p.gamma = gamma;
p.lambda = lambda; p.xi = xi; p.c = c; p.X0 = X0;
p.qK = qK; p.umin = umin; p.umax = umax;
p.maxIter = maxIter; p.tol = tol; p.omega = omega;
p.useClamping = useClamping;
p.dimY = n + 3*n*m;

%% Pre-entry phase: simulate only once, technologies 1 and 2 active
Y_pre = forwardSimulationPreEntry(y0,p);
yEntry = Y_pre(:,end);
[~,~,D_entry,~] = unpackTrajectory(Y_pre,n,m,T1);
gapT1 = mean(Y_pre(n+n+(1:n), end)) - mean(Y_pre(n+(1:n), end));
fprintf('At T1:  A2-A1 gap = %+.4f,  mean D1 = %.4f,  mean D2 = %.4f\n', ...
    gapT1, mean(D_entry(end,:,1)), mean(D_entry(end,:,2)));
fprintf('Target: gap>0 at T1, gap<0 at Tfin (unctrl), gap>0 at Tfin (ctrl)\n');

%% Baseline post-entry simulation with no control
Uzero_post = zeros(H,n,m);
Y_un_post = forwardSimulationPost(yEntry,Uzero_post,p);
Y_un = [Y_pre, Y_un_post(:,2:end)];
[S_un,A_un,D_un,X_un] = unpackTrajectory(Y_un,n,m,Tfin);

%% Problem 4: kingmaker control on technology 3 after entry
[Y_ctrl_post,U_post,MU_post,err,objVals,quadCostVals] = forwardBackwardSweepKingmakerFast(yEntry,p);
Y_ctrl = [Y_pre, Y_ctrl_post(:,2:end)];
U_ctrl = zeros(Tfin,n,m);
U_ctrl(T1+1:Tfin,:,:) = U_post;
[S_ctrl,A_ctrl,D_ctrl,X_ctrl] = unpackTrajectory(Y_ctrl,n,m,Tfin);

%% Diagnostics
fprintf('\n=== Problem 4: Kingmaker control, quadratic-cost version ===\n');
fprintf('Pre-entry gap at T1, uncontrolled: mean A2 - mean A1 = %.4f\n', ...
    mean(A_un(T1+1,:,2)) - mean(A_un(T1+1,:,1)));
fprintf('Final gap uncontrolled: mean A2 - mean A1 = %.4f\n', ...
    mean(A_un(Tfin+1,:,2)) - mean(A_un(Tfin+1,:,1)));
fprintf('Final gap controlled:   mean A2 - mean A1 = %.4f\n', ...
    mean(A_ctrl(Tfin+1,:,2)) - mean(A_ctrl(Tfin+1,:,1)));
fprintf('Quadratic control cost term = %.4f\n', computeQuadraticCostPost(U_post,p));
fprintf('Max control outside technology 3 = %.4e\n', max(max(max(abs(U_ctrl(:,:,[1,2]))))));
fprintf('Max control before entry = %.4e\n', max(max(max(abs(U_ctrl(1:T1,:,:))))));
fprintf('Max control at final ineffective step u(Tfin-1) = %.4e\n', max(max(abs(U_ctrl(Tfin,:,:)))));
printStateDiagnostics('uncontrolled',S_un,A_un,D_un,X_un);
printStateDiagnostics('controlled',S_ctrl,A_ctrl,D_ctrl,X_ctrl);
printPMPDiagnostics(U_post,p);

%% Plots
plotResultsProblem4(S_un,A_un,D_un,X_un,S_ctrl,A_ctrl,D_ctrl,X_ctrl,U_ctrl,err,quadCostVals,p);

%% ============================================================
%  Local functions
% ============================================================

function [Y,U,MU,err,objVals,quadCostVals] = forwardBackwardSweepKingmakerFast(yEntry,p)
    H = p.H; n = p.n; m = p.m;
    U = zeros(H,n,m);
    err = zeros(p.maxIter,1);
    objVals = zeros(p.maxIter,1);
    quadCostVals = zeros(p.maxIter,1);

    for h = 1:p.maxIter
        Y = forwardSimulationPost(yEntry,U,p);
        MU = backwardAdjointPost(Y,U,p);
        U_tilde = updateControlKingmakerPost(MU,p);

        U_new = p.omega*U_tilde + (1-p.omega)*U;
        U_new = enforceKingmakerAdmissibilityPost(U_new,p);

        err(h) = max(abs(U_new(:)-U(:)));
        objVals(h) = computeObjectiveKingmakerPost(Y,U,p);
        quadCostVals(h) = computeQuadraticCostPost(U_new,p);

        if mod(h,10)==0 || h==1
            fprintf('Iter %3d: err = %.3e, J = %.4f, quad cost = %.4f\n', ...
                h,err(h),objVals(h),quadCostVals(h));
        end

        U = U_new;
        if err(h) < p.tol
            err = err(1:h); objVals = objVals(1:h); quadCostVals = quadCostVals(1:h);
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
    % Quadratic-cost PMP maximizer: max b*u - qK*u^2 over [umin,umax]
    % => u = proj_{[umin,umax]}( b / (2*qK) ), b = c*(1-X0)*mu_X(t+1).
    H = p.H; n = p.n; m = p.m; k = p.k_new;
    U_tilde = zeros(H,n,m);

    coeff  = p.c(:,k) .* (1 - p.X0(:,k));          % n x 1
    idxXk  = n + 2*n*m + (k-1)*n + (1:n);          % X^[k] rows in state vector
    b      = (MU(idxXk, 2:H) .* coeff)';           % (H-1) x n

    U_tilde(1:H-1, :, k) = min(p.umax, max(p.umin, b / (2*p.qK)));

    U_tilde = enforceKingmakerAdmissibilityPost(U_tilde,p);
end

function U = enforceKingmakerAdmissibilityPost(U,p)
    % Post-entry U has rows corresponding to global times T1,...,Tfin-1.
    % Control only technology 3, and not at final ineffective step.
    U(:,:,1) = 0;
    U(:,:,2) = 0;
    U(p.H,:,:) = 0;      % global math t = Tfin-1
    U = min(max(U,p.umin),p.umax);
end

function J = computeObjectiveKingmakerPost(Y,U,p)
    n = p.n;
    gap = sum(Y(n + n + (1:n), p.H+1)) - sum(Y(n + (1:n), p.H+1));
    J   = gap - computeQuadraticCostPost(U,p);
end

function C = computeQuadraticCostPost(U,p)
    C = p.qK * sum(U(1:p.H-1, :, p.k_new).^2, 'all');
end

function Bused = computeLinearBudgetPost(U,p)
    % Diagnostic: total linear effort over global t=T1,...,Tfin-2.
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

function plotResultsProblem4(S_un,A_un,D_un,X_un,S_ctrl,A_ctrl,D_ctrl,X_ctrl,U_ctrl,err,quadCostVals,p) %#ok<INUSD>
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
    xline(T1,':k','Interpreter','latex','LineWidth',1.3);
    xlabel('Time','Interpreter','latex'); ylabel('Mean adoption','Interpreter','latex');
    %title('Problem 4: adoption trajectories','Interpreter','latex');
    grid on; set(gca,'FontSize',14);
    legend({'Unctrl $A^{[1]}$','Ctrl $A^{[1]}$','Unctrl $A^{[2]}$','Ctrl $A^{[2]}$', ...
            'Unctrl $A^{[3]}$','Ctrl $A^{[3]}$'},'Interpreter','latex','Location','north','NumColumns',3);

    nexttile; hold on;
    gap_un   = meanA_un(:,2)   - meanA_un(:,1);
    gap_ctrl = meanA_ctrl(:,2) - meanA_ctrl(:,1);
    fill([timeState, fliplr(timeState)], [gap_un', fliplr(gap_ctrl')], ...
        [0.2 0.6 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none');   % shaded improvement
    plot(timeState, gap_un,   '--k', 'LineWidth', 1.8);
    plot(timeState, gap_ctrl, '-',   'LineWidth', 2.4, 'Color', [0 0.45 0.74]);
    yline(0,   ':', 'LineWidth', 1.2);
    xline(T1,  ':', 'LineWidth', 1.3, 'Color', 'k');
    xlabel('Time','Interpreter','latex');
    ylabel('$\overline A^{[2]}-\overline A^{[1]}$','Interpreter','latex');
    %title('Aggregate rank gap','Interpreter','latex');
    legend('Kingmaker gain','Uncontrolled','Controlled','Interpreter','latex','Location','best');
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
    %title('Dissatisfied populations','Interpreter','latex');
    grid on; set(gca,'FontSize',14);

    nexttile; hold on;
    plot(timeCtrl, meanU3, 'Color', colorsTech(3,:), 'LineWidth', 2.2);
    xline(T1, ':', 'LineWidth', 1.3, 'Color', 'k');
    xlabel('Time','Interpreter','latex');
    ylabel('Mean $u^{[3]}$','Interpreter','latex');
    title('Kingmaker control profile','Interpreter','latex');
    grid on; set(gca,'FontSize',14);

    % Extra figure: gap difference (controlled - uncontrolled) — makes even
    % small kingmaker gains clearly visible.
    figure; hold on;
    delta_gap = gap_ctrl - gap_un;
    area(timeState, delta_gap, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(timeState, delta_gap, 'Color', [0 0.3 0.7], 'LineWidth', 2.2);
    yline(0, '-k', 'LineWidth', 1);
    xline(T1, ':', 'LineWidth', 1.3, 'Color', 'k', 'Label', 'Entry of tech 3', ...
        'Interpreter', 'latex');
    xlabel('Time','Interpreter','latex');
    ylabel('$\Delta\mathrm{gap} = \overline A^{[2]}_{\mathrm{ctrl}} - \overline A^{[2]}_{\mathrm{unctrl}}$','Interpreter','latex');
    title('Kingmaker effect on rank gap (controlled $-$ uncontrolled)','Interpreter','latex');
    grid on; set(gca,'FontSize',14);

    figure;
    tiledlayout(2,1);
    nexttile;
    semilogy(err,'LineWidth',2);
    xlabel('Iteration','Interpreter','latex'); ylabel('Control update error','Interpreter','latex');
    title('Forward-backward sweep convergence','Interpreter','latex');
    grid on; set(gca,'FontSize',14);

    nexttile;
    plot(quadCostVals,'LineWidth',2);
    xlabel('Iteration','Interpreter','latex'); ylabel('Quadratic control cost','Interpreter','latex');
    title('Quadratic intervention cost over iterations','Interpreter','latex');
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
    fprintf('PMP check: linear effort = %.4f, quadratic cost = %.4f, fractional active controls = %d\n', ...
        computeLinearBudgetPost(U,p),computeQuadraticCostPost(U,p),numel(frac));
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
