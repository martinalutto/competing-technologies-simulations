%% run_control_three_problems_noD.m
% Two optimal-control problems for the adoption-opinion model.
%
% Problem 1:
%   Cumulative system-level adoption maximization.
%   Control acts on all technologies and is penalized in the objective.
%
% Problem 2:
%   Targeted terminal adoption maximization under a budget constraint.
%   Control acts only on the target technology k0 and must satisfy
%       sum_t u(t)' Q u(t) <= B.
%
% Time convention:
%   State arrays S,A,D,X have size (T+1) x n x m.
%   Row 1 corresponds to mathematical time t = 0.
%   Control array U has size T x n x m.
%   Row 1 corresponds to mathematical time t = 0.

clear; close all; clc;
rng(2);

%% ============================================================
%  PARAMETERS
% ============================================================

n = 20;              % number of agents
m = 5;               % number of technologies/layers
T = 300;             % time horizon
k0 = 1;              % target technology for Problems 2 and 3

%% Model parameters
beta  = zeros(n,m);
delta = zeros(n,m);
gamma = zeros(n,m);

beta(:,1)  = 0.05 + 0.05*rand(n,1);
beta(:,2)  = 0.10 + 0.05*rand(n,1);
beta(:,3)  = 0.18 + 0.01*rand(n,1);
beta(:,4)  = 0.10 + 0.05*rand(n,1);
beta(:,5)  = 0.05 + 0.03*rand(n,1);

delta(:,1) = 0.08 + 0.04*rand(n,1);
delta(:,2) = 0.10 + 0.04*rand(n,1);
delta(:,3) = 0.15 + 0.05*rand(n,1);
delta(:,4) = 0.20 + 0.08*rand(n,1);
delta(:,5) = 0.18 + 0.10*rand(n,1);

gamma(:,1) = 0.10 + 0.04*rand(n,1);
gamma(:,2) = 0.12 + 0.04*rand(n,1);
gamma(:,3) = 0.14 + 0.04*rand(n,1);
gamma(:,4) = 0.12 + 0.04*rand(n,1);
gamma(:,5) = 0.14 + 0.04*rand(n,1);

lambda = 0.30*ones(n,m);
xi     = 0.10*ones(n,m);
c      = 1 - lambda - xi;

%% Network matrices
W = rand(n,n);
W = W ./ sum(W,2);

tilde_W = rand(n,n);
tilde_W = tilde_W ./ sum(tilde_W,2);

%% Fixed initial opinions used in x_U
% x_i^{[k]}(0) in (0,1] (Assumption): strictly positive lower bound.
X0 = 0.05 + 0.65*rand(n,m);

%% Objective parameters

% Problem 1: weights on cumulative adoption of all technologies
r1 = [0.6, 0.7, 0.50, 0.4, 0.50];

% Problem 2: terminal target adoption reward. The objective is
% J_2 = 1' a^[k0](T), hence r2 = 1 (terminal costate on a^[k0] equals 1).
% Dissatisfied variables are not included in the objective; they are only
% monitored in the plots as a side effect of the intervention.
r2 = 1.0;
r3 = 1.0;   % unused (kept only to preserve function signatures)

% Budget for Problem 2.
% If umax = 1 and B is very small, the target-budget control naturally
% tends to concentrate near the end. Increase B to obtain a less impulsive
% campaign, or reduce umax.
B  = 40;

% Control cost matrix Q = diag(q_i^[k]) represented as an n-by-m matrix
q = 0.08*ones(n,m);

umin = 0;
umax = 1;

%% Forward-backward sweep parameters
maxIter = 200;
tol = 1e-4;
omega = 0.30;

%% ============================================================
%  INITIAL CONDITIONS
% ============================================================

S0 = zeros(1,n);
A0 = zeros(1,n,m);
D0 = zeros(1,n,m);
Xinit = zeros(1,n,m);

% Example initialization.
% You can modify these values depending on the scenario.
A0(1,:,1) = 0.05*rand(1,n);
A0(1,:,2) = 0.03*rand(1,n);
A0(1,:,3) = 0.02*rand(1,n);
A0(1,:,4) = 0.01*rand(1,n);
A0(1,:,5) = 0.01*rand(1,n);

D0(1,:,:) = 0;
S0(1,:) = max(0, 1 - squeeze(sum(A0,3))' - squeeze(sum(D0,3))');

% Initial opinion state equals the fixed initial opinions.
for k = 1:m
    Xinit(1,:,k) = X0(:,k)';
end

%% ============================================================
%  UNCONTROLLED SIMULATION
% ============================================================

Uzero = zeros(T,n,m);

[S_un,A_un,D_un,X_un] = forwardSimulation( ...
    S0,A0,D0,Xinit,Uzero, ...
    W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,T,n,m);

%% ============================================================
%  PROBLEM 1: CUMULATIVE SYSTEM-LEVEL CONTROL
% ============================================================

[S_cum,A_cum,D_cum,X_cum,U_cum,err_cum,eta_cum] = forwardBackwardSweep( ...
    S0,A0,D0,Xinit, ...
    W,tilde_W,beta,delta,gamma,lambda,xi,c,X0, ...
    r1,r2,r3,q,B,k0,umin,umax,T,n,m,maxIter,tol,omega, ...
    "cumulative");

%% ============================================================
%  PROBLEM 2: TARGETED TERMINAL CONTROL WITH BUDGET
% ============================================================

[S_tar,A_tar,D_tar,X_tar,U_tar,err_tar,eta_tar] = forwardBackwardSweep( ...
    S0,A0,D0,Xinit, ...
    W,tilde_W,beta,delta,gamma,lambda,xi,c,X0, ...
    r1,r2,r3,q,B,k0,umin,umax,T,n,m,maxIter,tol,omega, ...
    "target_budget");

%% ============================================================
%  OBJECTIVES AND DIAGNOSTICS
% ============================================================

J1_un  = computeObjectiveProblem1(A_un,Uzero,r1,q,T,n,m);
J1_cum = computeObjectiveProblem1(A_cum,U_cum,r1,q,T,n,m);

J2_un  = computeObjectiveProblem2(A_un,D_un,r2,r3,k0,T);
J2_tar = computeObjectiveProblem2(A_tar,D_tar,r2,r3,k0,T);

budget_tar = computeBudget(U_tar,q,T,n,m);

fprintf('\n================ RESULTS ================\n');
fprintf('Problem 1 objective, uncontrolled: %.6f\n', J1_un);
fprintf('Problem 1 objective, controlled:   %.6f\n', J1_cum);
fprintf('Problem 2 objective, uncontrolled: %.6f\n', J2_un);
fprintf('Problem 2 objective, controlled:   %.6f\n', J2_tar);
fprintf('Problem 2 budget used:             %.6f / %.6f\n', budget_tar, B);
fprintf('Problem 2 eta estimate:            %.6e\n', eta_tar);

U_other = U_tar;
U_other(:,:,k0) = 0;
fprintf('Max control outside target layer in Problem 2: %.3e\n', max(abs(U_other(:))));

checkStateBounds(S_un,A_un,D_un,X_un,"uncontrolled");
checkStateBounds(S_cum,A_cum,D_cum,X_cum,"problem 1");
checkStateBounds(S_tar,A_tar,D_tar,X_tar,"problem 2");

%% ============================================================
%  PLOTS
% ============================================================

timeState = 0:size(A_un,1)-1;
timeCtrl_cum = 0:size(U_cum,1)-1;
timeCtrl_tar = 0:size(U_tar,1)-1;

%% Medie sui nodi
meanA_un  = squeeze(mean(A_un,  2));   % (T+1) x m
meanD_un  = squeeze(mean(D_un,  2));   % (T+1) x m

meanA_cum = squeeze(mean(A_cum, 2));   % (T+1) x m
meanD_cum = squeeze(mean(D_cum, 2));   % (T+1) x m

meanA_tar = squeeze(mean(A_tar, 2));     % (T+1) x m
meanD_tar = squeeze(mean(D_tar, 2));     % (T+1) x m

%% Palette: same color for A and D of the same technology
colorsTech = [
    0.00 0.45 0.74   % technology 1
    0.85 0.33 0.10   % technology 2
    0.47 0.67 0.19   % technology 3
    0.49 0.18 0.56   % technology 4
    0.93 0.69 0.13   % technology 5
];
if m > size(colorsTech,1)
    colorsTech = lines(m);
end
colorsA = colorsTech;
colorsD = colorsTech;

%% ============================================================
%  PROBLEM 1: all technologies, controlled vs uncontrolled
% ============================================================

figure;
tiledlayout(2,1);

% --- Adoption, all technologies ---
nexttile; hold on;

for k = 1:m
    plot(timeState, meanA_un(:,k), '--', 'LineWidth', 1.6, 'Color', colorsA(k,:));
    plot(timeState, meanA_cum(:,k), '-', 'LineWidth', 2.2, 'Color', colorsA(k,:));
end

xlabel('Time', 'Interpreter','latex');
ylabel('Mean adoption', 'Interpreter','latex');
title('Problem 1: adoption of all technologies', 'Interpreter','latex');
grid on;
set(gca,'FontSize',14);

legA = cell(1,2*m);
for k = 1:m
    legA{2*k-1} = ['Unctrl $A^{[',num2str(k),']}$'];
    legA{2*k}   = ['Ctrl $A^{[',num2str(k),']}$'];
end
legend(legA, 'Interpreter','latex', 'Location','bestoutside', 'NumColumns',2);

% --- Dissatisfied, all technologies ---
nexttile; hold on;

for k = 1:m
    plot(timeState, meanD_un(:,k), '--', 'LineWidth', 1.6, 'Color', colorsD(k,:));
    plot(timeState, meanD_cum(:,k), '-', 'LineWidth', 2.2, 'Color', colorsD(k,:));
end

xlabel('Time', 'Interpreter','latex');
ylabel('Mean dissatisfied', 'Interpreter','latex');
title('Problem 1: dissatisfied population for all technologies', 'Interpreter','latex');
grid on;
set(gca,'FontSize',14);

legD = cell(1,2*m);
for k = 1:m
    legD{2*k-1} = ['Unctrl $D^{[',num2str(k),']}$'];
    legD{2*k}   = ['Ctrl $D^{[',num2str(k),']}$'];
end
legend(legD, 'Interpreter','latex', 'Location','bestoutside', 'NumColumns',2);

%% ============================================================
%  PROBLEM 2: target technology only, controlled vs uncontrolled
% ============================================================

figure;
tiledlayout(2,1);

% --- Target adoption ---
nexttile; hold on;

plot(timeState, meanA_un(:,k0), '--k', 'LineWidth', 1.8);
plot(timeState, meanA_tar(:,k0), '-',  'LineWidth', 2.4, ...
    'Color', colorsA(k0,:));

xlabel('Time', 'Interpreter','latex');
ylabel('Mean adoption', 'Interpreter','latex');
title(['Problem 2: target adoption, $k_0=',num2str(k0),'$'], ...
    'Interpreter','latex');
legend('Uncontrolled', 'Controlled', ...
    'Interpreter','latex', ...
    'Location','best');
grid on;
set(gca,'FontSize',14);

% --- Target dissatisfied ---
nexttile; hold on;

plot(timeState, meanD_un(:,k0), '--k', 'LineWidth', 1.8);
plot(timeState, meanD_tar(:,k0), '-',  'LineWidth', 2.4, ...
    'Color', colorsD(k0,:));

xlabel('Time', 'Interpreter','latex');
ylabel('Mean dissatisfied', 'Interpreter','latex');
title(['Problem 2: target dissatisfied, $k_0=',num2str(k0),'$'], ...
    'Interpreter','latex');
legend('Uncontrolled', 'Controlled', ...
    'Interpreter','latex', ...
    'Location','best');
grid on;
set(gca,'FontSize',14);

%% ============================================================
%  CONTROLS OVER TIME
% ============================================================

% Problem 1: average control over all nodes and technologies
meanU_cum_all = squeeze(mean(mean(U_cum,2),3));  % T x 1

% Problem 1: average control by technology
meanU_cum_byTech = squeeze(mean(U_cum,2));       % T x m

% Problem 2: average control on the target technology only
meanU_tar_target = squeeze(mean(U_tar(:,:,k0),2)); % T x 1

figure;
tiledlayout(2,1);

% --- Overall comparison between the two controls ---
nexttile; hold on;

plot(timeCtrl_cum, meanU_cum_all, 'LineWidth', 2.2);
plot(timeCtrl_tar, meanU_tar_target, 'LineWidth', 2.2);

xlabel('Time', 'Interpreter','latex');
ylabel('Mean control', 'Interpreter','latex');
title('Average control over time', 'Interpreter','latex');
legend('Problem 1: mean over all technologies', ...
       'Problem 2: mean on target technology', ...
       'Interpreter','latex', ...
       'Location','best');
grid on;
set(gca,'FontSize',14);

% --- Problem 1: control allocation across technologies ---
nexttile; hold on;

for k = 1:m
    plot(timeCtrl_cum, meanU_cum_byTech(:,k), ...
        'LineWidth', 2.0);
end

xlabel('Time', 'Interpreter','latex');
ylabel('Mean control by technology', 'Interpreter','latex');
title('Problem 1: control allocation across technologies', ...
    'Interpreter','latex');

legend_entries = cell(1,m);
for k = 1:m
    legend_entries{k} = ['$u^{[',num2str(k),']}$'];
end
legend(legend_entries, ...
    'Interpreter','latex', ...
    'Location','best', ...
    'NumColumns',m);
grid on;
set(gca,'FontSize',14);

%% ============================================================
%  CONVERGENCE
% ============================================================

figure; hold on;
semilogy(err_cum, 'LineWidth', 2);
semilogy(err_tar, 'LineWidth', 2);
xlabel('Iteration', 'Interpreter','latex');
ylabel('$\max_t \|u^{(h+1)}(t)-u^{(h)}(t)\|$', 'Interpreter','latex');
legend('Problem 1', 'Problem 2', 'Interpreter','latex', 'Location','best');
title('Forward--backward sweep convergence', 'Interpreter','latex');
grid on;
set(gca,'FontSize',14);

%% ============================================================
%  FUNCTIONS
% ============================================================

function [S,A,D,X,U,err,etaLast] = forwardBackwardSweep( ...
    S0,A0,D0,Xinit, ...
    W,tilde_W,beta,delta,gamma,lambda,xi,c,X0, ...
    r1,r2,r3,q,B,k0,umin,umax,T,n,m,maxIter,tol,omega, ...
    problemType)

    U = zeros(T,n,m);
    err = zeros(maxIter,1);
    etaLast = NaN;

    for h = 1:maxIter

        [S,A,D,X] = forwardSimulation( ...
            S0,A0,D0,Xinit,U, ...
            W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,T,n,m);

        MU = backwardAdjoint( ...
            S,A,D,X,U, ...
            W,tilde_W,beta,delta,gamma,lambda,xi,c,X0, ...
            r1,r2,r3,q,k0,T,n,m,problemType);

        switch problemType
            case "cumulative"
                U_tilde = updateControlCumulative( ...
                    MU,c,X0,q,umin,umax,T,n,m);
                etaLast = NaN;

            case "target_budget"
                [U_tilde,etaLast] = updateControlTargetBudget( ...
                    MU,c,X0,q,B,k0,umin,umax,T,n,m);

            otherwise
                error('Unknown problem type.');
        end

        U_new = omega*U_tilde + (1-omega)*U;

        err(h) = maxControlDiff(U_new,U);

        U = U_new;

        if mod(h,10)==0 || h==1
            fprintf('%s iter %4d: err = %.4e\n', problemType, h, err(h));
        end

        if err(h) < tol
            err = err(1:h);
            fprintf('%s converged at iter %d with err %.4e\n', problemType, h, err(end));
            break;
        end

        if h == maxIter
            fprintf('%s reached maxIter with err %.4e\n', problemType, err(h));
        end
    end

    [S,A,D,X] = forwardSimulation( ...
        S0,A0,D0,Xinit,U, ...
        W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,T,n,m);
end

function [S,A,D,X] = forwardSimulation( ...
    S0,A0,D0,Xinit,U, ...
    W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,T,n,m)

    S = zeros(T+1,n);
    A = zeros(T+1,n,m);
    D = zeros(T+1,n,m);
    X = zeros(T+1,n,m);

    S(1,:) = S0;
    A(1,:,:) = A0;
    D(1,:,:) = D0;
    X(1,:,:) = Xinit;

    for t = 1:T
        Sold = S(t,:)';
        Aold = squeeze(A(t,:,:));
        Dold = squeeze(D(t,:,:));
        Xold = squeeze(X(t,:,:));
        Uold = squeeze(U(t,:,:));

        [Snew,Anew,Dnew,Xnew] = oneStepState( ...
            Sold,Aold,Dold,Xold,Uold, ...
            W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,n,m);

        S(t+1,:) = Snew';
        A(t+1,:,:) = Anew;
        D(t+1,:,:) = Dnew;
        X(t+1,:,:) = Xnew;
    end
end

function [Snew,Anew,Dnew,Xnew] = oneStepState( ...
    S,A,D,X,U, ...
    W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,n,m)

    Wa = W * A;                                   % n x m  (matrix product)
    totalForce  = sum(beta .* X .* Wa, 2);
    Snew        = S .* (1 - totalForce);

    Anew = zeros(n,m);
    Dnew = zeros(n,m);
    Xnew = zeros(n,m);

    totalGammaX = sum(gamma.*X, 2);
    totalD      = sum(D, 2);

    for k = 1:m
        crossD       = totalD      - D(:,k);
        switchingOut = totalGammaX - gamma(:,k).*X(:,k);

        Anew(:,k) = A(:,k) + beta(:,k).*X(:,k).*S.*Wa(:,k) ...
            - delta(:,k).*A(:,k) + gamma(:,k).*X(:,k).*crossD;

        Dnew(:,k) = D(:,k) .* (1 - switchingOut) + delta(:,k).*A(:,k);

        XU         = X0(:,k).*(1 - U(:,k)) + U(:,k);
        Xnew(:,k)  = c(:,k).*XU + lambda(:,k).*(tilde_W*X(:,k)) + xi(:,k).*Wa(:,k);
    end
end

function MU = backwardAdjoint( ...
    S,A,D,X,U, ...
    W,tilde_W,beta,delta,gamma,lambda,xi,c,X0, ...
    r1,r2,r3,q,k0,T,n,m,problemType)

    dimY = n + 3*n*m;
    MU = zeros(dimY,T+1);

    % Terminal condition
    if problemType == "target_budget"
        % Phi_2 = 1' a^[k0](T)  ->  terminal costate on a^[k0] equals r2 (=1)
        MU(n + (k0-1)*n + (1:n), T+1) = r2;
    elseif problemType ~= "cumulative"
        error('Unknown problem type: %s', problemType);
    end
    % cumulative: no terminal reward, MU(:,T+1) stays zero

    gradL = runningGradientY(r1,k0,n,m,problemType);   % constant over t

    for t = T:-1:1

        y = packState(S(t,:)', squeeze(A(t,:,:)), squeeze(D(t,:,:)), squeeze(X(t,:,:)), n,m);
        u = squeeze(U(t,:,:));

        Jy = analyticJacobianY( ...
            y,u,W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,n,m);

        MU(:,t) = gradL + Jy'*MU(:,t+1);
    end
end

function gradL = runningGradientY(r1,k0,n,m,problemType)
    dimY  = n + 3*n*m;
    gradL = zeros(dimY,1);
    if problemType == "cumulative"
        for k = 1:m
            gradL(n + (k-1)*n + (1:n)) = r1(k);
        end
    end
    % target_budget: no running reward, gradL stays zero
end

function U_tilde = updateControlCumulative(MU,c,X0,q,umin,umax,T,n,m)
    U_tilde = zeros(T,n,m);
    for k = 1:m
        coeff  = c(:,k) .* (1 - X0(:,k));           % n x 1
        idxXk  = n + 2*n*m + (k-1)*n + (1:n);
        b      = (MU(idxXk, 2:T+1) .* coeff)';      % T x n
        uCand  = b ./ (2 * q(:,k)');                 % T x n  (broadcasting)
        U_tilde(:,:,k) = min(umax, max(umin, uCand));
    end
end

function [U_tilde,eta] = updateControlTargetBudget( ...
    MU,c,X0,q,B,k0,umin,umax,T,n,m)

    % Control is allowed only on layer k0.
    % We search eta so that sum_t u'Q u <= B.
    % For eta > 0:
    %   u_i(t) = Proj_[umin,umax] b_i(t)/(2 eta q_i)

    coeff  = c(:,k0) .* (1 - X0(:,k0));
    idxXk0 = n + 2*n*m + (k0-1)*n + (1:n);
    b      = (MU(idxXk0, 2:T+1) .* coeff)';   % T x n

    % If all marginal benefits are nonpositive, the optimal control is zero.
    if max(b(:)) <= 0
        U_tilde = zeros(T,n,m);
        eta = Inf;
        return;
    end

    function Utmp = controlForEta(etaValue)
        Utmp = zeros(T,n,m);
        uCand = b ./ (2*etaValue * q(:,k0)');   % T x n  (broadcasting)
        Utmp(:,:,k0) = min(umax, max(umin, uCand));
    end

    % Very small eta approximates the bound-constrained linear maximizer.
    etaLow = 1e-12;
    U_low = controlForEta(etaLow);
    budgetLow = computeBudget(U_low,q,T,n,m);

    % If even the saturated control respects the budget, then the budget is
    % inactive in practice. We return the saturated pointwise maximizer.
    if budgetLow <= B
        U_tilde = U_low;
        eta = etaLow;
        return;
    end

    % Otherwise find etaHigh such that budget <= B.
    etaHigh = 1;
    U_high = controlForEta(etaHigh);
    budgetHigh = computeBudget(U_high,q,T,n,m);

    while budgetHigh > B
        etaHigh = 2*etaHigh;
        U_high = controlForEta(etaHigh);
        budgetHigh = computeBudget(U_high,q,T,n,m);

        if etaHigh > 1e12
            warning('Could not bracket eta for budget search.');
            break;
        end
    end

    % Bisection on eta.
    for iter = 1:70
        etaMid = 0.5*(etaLow + etaHigh);
        U_mid = controlForEta(etaMid);
        budgetMid = computeBudget(U_mid,q,T,n,m);

        if budgetMid > B
            etaLow = etaMid;
        else
            etaHigh = etaMid;
        end
    end

    eta = etaHigh;
    U_tilde = controlForEta(eta);
end

function J = analyticJacobianY( ...
    y,u,W,tilde_W,beta,delta,gamma,lambda,xi,c,X0,n,m)

    dimY = n + 3*n*m;
    J = zeros(dimY,dimY);

    [S,A,D,X] = unpackState(y,n,m);

    Wa = zeros(n,m);
    for k = 1:m
        Wa(:,k) = W*A(:,k);
    end

    totalForce = zeros(n,1);
    for k = 1:m
        totalForce = totalForce + beta(:,k).*X(:,k).*Wa(:,k);
    end

    totalGammaX = sum(gamma.*X,2);
    totalD = sum(D,2);

    %% Derivatives of S_next
    for i = 1:n
        rowS = indexS(i);

        % d S_next_i / d S_i
        J(rowS,indexS(i)) = 1 - totalForce(i);

        for k = 1:m
            % d S_next_i / d X_i^k
            J(rowS,indexX(i,k,n,m)) = J(rowS,indexX(i,k,n,m)) ...
                - S(i)*beta(i,k)*Wa(i,k);

            % d S_next_i / d A_l^k through Wa_i^k
            for l = 1:n
                J(rowS,indexA(l,k,n)) = J(rowS,indexA(l,k,n)) ...
                    - S(i)*beta(i,k)*X(i,k)*W(i,l);
            end
        end
    end

    %% Derivatives of A_next, D_next, X_next
    for i = 1:n
        for k = 1:m

            rowA = indexA(i,k,n);
            rowD = indexD(i,k,n,m);
            rowX = indexX(i,k,n,m);

            crossD = totalD(i) - D(i,k);
            switchingOut = totalGammaX(i) - gamma(i,k)*X(i,k);

            %% A_next_i^k
            % wrt S_i
            J(rowA,indexS(i)) = beta(i,k)*X(i,k)*Wa(i,k);

            % wrt A_l^k through Wa_i^k
            for l = 1:n
                J(rowA,indexA(l,k,n)) = J(rowA,indexA(l,k,n)) ...
                    + beta(i,k)*X(i,k)*S(i)*W(i,l);
            end

            % wrt A_i^k from A_i^k - delta A_i^k
            J(rowA,indexA(i,k,n)) = J(rowA,indexA(i,k,n)) + 1 - delta(i,k);

            % wrt D_i^h, h ~= k
            for h = 1:m
                if h ~= k
                    J(rowA,indexD(i,h,n,m)) = J(rowA,indexD(i,h,n,m)) ...
                        + gamma(i,k)*X(i,k);
                end
            end

            % wrt X_i^k
            J(rowA,indexX(i,k,n,m)) = beta(i,k)*S(i)*Wa(i,k) ...
                + gamma(i,k)*crossD;

            %% D_next_i^k
            % wrt D_i^k
            J(rowD,indexD(i,k,n,m)) = 1 - switchingOut;

            % wrt A_i^k
            J(rowD,indexA(i,k,n)) = delta(i,k);

            % wrt X_i^h, h ~= k
            for h = 1:m
                if h ~= k
                    J(rowD,indexX(i,h,n,m)) = J(rowD,indexX(i,h,n,m)) ...
                        - D(i,k)*gamma(i,h);
                end
            end

            %% X_next_i^k
            % wrt A_l^k
            for l = 1:n
                J(rowX,indexA(l,k,n)) = J(rowX,indexA(l,k,n)) ...
                    + xi(i,k)*W(i,l);
            end

            % wrt X_l^k
            for l = 1:n
                J(rowX,indexX(l,k,n,m)) = J(rowX,indexX(l,k,n,m)) ...
                    + lambda(i,k)*tilde_W(i,l);
            end

            % No derivative wrt U here: this is Jacobian wrt state.
            % The control dependence is used explicitly in updateControl.
        end
    end
end

function y = packState(S,A,D,X,n,m)
    y = [S(:); A(:); D(:); X(:)];
end

function [S,A,D,X] = unpackState(y,n,m)

    S = y(1:n);

    startA = n + 1;
    endA = n + n*m;
    A = reshape(y(startA:endA),n,m);

    startD = endA + 1;
    endD = endA + n*m;
    D = reshape(y(startD:endD),n,m);

    startX = endD + 1;
    endX = endD + n*m;
    X = reshape(y(startX:endX),n,m);
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

function val = maxControlDiff(Unew,Uold)
    d   = reshape(Unew - Uold, size(Unew,1), []);   % T x (n*m)
    val = max(sqrt(sum(d.^2, 2)));
end

function Bused = computeBudget(U,q,T,n,m) %#ok<INUSD>
    Bused = sum(q .* squeeze(sum(U.^2, 1)), 'all');
end

function J = computeObjectiveProblem1(A,U,r1,q,T,n,m) %#ok<INUSD>
    r1mat          = reshape(r1, 1, 1, []);                  % 1 x 1 x m
    adoptionReward = sum(r1mat .* A(1:T,:,:), 'all');
    controlCost    = sum(q .* squeeze(sum(U.^2, 1)), 'all');
    J              = adoptionReward - controlCost;
end

function J = computeObjectiveProblem2(A,D,r2,r3,k0,T)
    %#ok<INUSD> D and r3 are intentionally unused: dissatisfied variables are monitored but not optimized.
    % Terminal state is row T+1.  J_2 = r2 * 1' a^[k0](T), with r2 = 1.
    J = r2*sum(A(T+1,:,k0));
end

function checkStateBounds(S,A,D,X,name)
    minS = min(S(:)); maxS = max(S(:));
    minA = min(A(:)); maxA = max(A(:));
    minD = min(D(:)); maxD = max(D(:));
    minX = min(X(:)); maxX = max(X(:));

    massErr = zeros(size(S,1),1);
    for t = 1:size(S,1)
        At = squeeze(A(t,:,:));
        Dt = squeeze(D(t,:,:));
        massErr(t) = max(abs(S(t,:)' + sum(At,2) + sum(Dt,2) - 1));
    end

    fprintf('\nState bounds: %s\n', name);
    fprintf('  S in [%.4f, %.4f]\n', minS, maxS);
    fprintf('  A in [%.4f, %.4f]\n', minA, maxA);
    fprintf('  D in [%.4f, %.4f]\n', minD, maxD);
    fprintf('  X in [%.4f, %.4f]\n', minX, maxX);
    fprintf('  max mass error = %.3e\n', max(massErr));
end
