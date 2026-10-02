% Original exploratory script: MATE dynamics with m = 5 technologies
% (per-agent loop; see src/simulate_adoption_multi.m for the vectorized version).
%% Parameters
rng(3);
n = 30;          % agents
m = 5;           % technologies
T = 200;         % horizon

%% Technology parameters (n x m)
beta   = rand(n,m);
beta   = beta ./ max(1, sum(beta,2));   % Assumption 1(iii)
delta  = rand(n,m);
gamma  = 0.5*rand(n,m);
gamma  = gamma ./ max(1, sum(gamma,2)); % Assumption 1(iii)

lambda = rand(n,m);
lambda = 0.8*lambda ./ sum(lambda,2);   % rows sum to 0.8
xi = rand(n,m) .* (1 - lambda);

%% Networks
W = rand(n,n);
W = W ./ sum(W,2);

tilde_W = rand(n,n);
tilde_W = tilde_W ./ sum(tilde_W,2);

S = zeros(T,n);
A = zeros(T,n,m);
D = zeros(T,n,m);
X = zeros(T,n,m);

%% Initialization
A(1,:,1) = rand(1,n);
A(1,:,2) = 0.8*rand(1,n);
A(1,:,3) = zeros(1,n);
A(1,:,4) = zeros(1,n);
A(1,:,5) = 0.1*rand(1,n);
S(1,:) = max(zeros(1,n),ones(1,n)-sum(A(1,:,:),3));
D(1,:,:) = 0;
X(1,:,:) = 0.4;

%% Simulation
for t = 1:T-1
    
    % Wa(i,k) = sum_j W(i,j) * A(t,j,k)
    Wa = zeros(n,m);
    for k = 1:m
        Wa(:,k) = W * squeeze(A(t,:,k))';
    end
    
    for i = 1:n
        s_old = S(t,i);
        adoption_force = 0;
        for k = 1:m
            % Opinions
            X0ik = X(1,i,k);
            X_neighbors = tilde_W(i,:) * squeeze(X(t,:,k))';
            X(t+1,i,k) = (1 - lambda(i,k) - xi(i,k)) * X0ik + lambda(i,k) * X_neighbors + xi(i,k) * Wa(i,k);
            X(t+1,i,k) = max(0, min(1, X(t+1,i,k)));

            % Susceptibles
            adoption_force = adoption_force + beta(i,k) * X(t,i,k) * Wa(i,k);
            S(t+1,i) = s_old - s_old * adoption_force;
            S(t+1,i) = max(0, min(1, S(t+1,i)));
        
            % Adopters
            dissatisfied_others = 0;
            for h = 1:m
                if h ~= k
                    dissatisfied_others = dissatisfied_others + D(t,i,h);
                end
            end

            A(t+1,i,k) = A(t,i,k) + beta(i,k) * X(t,i,k) * s_old * Wa(i,k) - delta(i,k) * A(t,i,k) + gamma(i,k) * X(t,i,k) * dissatisfied_others;
            A(t+1,i,k) = max(0, min(1, A(t+1,i,k)));

            % Dissatisfied
            switching_out = 0;
            for h = 1:m
                if h ~= k
                    switching_out = switching_out + gamma(i,h) * X(t,i,h);
                end
            end
          
            D(t+1,i,k) = D(t,i,k) + delta(i,k) * A(t,i,k) - D(t,i,k) * switching_out;
            D(t+1,i,k) = max(0, min(1, D(t+1,i,k)));
        end
    end
end

%% Node averages
meanS = mean(S,2);

meanA = zeros(T,m);
meanD = zeros(T,m);
meanX = zeros(T,m);

for k = 1:m
    meanA(:,k) = mean(squeeze(A(:,:,k)), 2);
    meanD(:,k) = mean(squeeze(D(:,:,k)), 2);
    meanX(:,k) = mean(squeeze(X(:,:,k)), 2);
end

[meanA(end,1), meanA(end,2),meanA(end,3)]

%% Plot adopters (warm colors) and dissatisfied (cool colors)
figure; hold on;
colorsA = [
    0.75 0.33 0.25
    0.86 0.49 0.20
    0.93 0.69 0.13
    0.80 0.52 0.25
    0.91 0.57 0.33
    0.70 0.42 0.18
];
colorsD = [
    0.16 0.44 0.52
    0.24 0.60 0.56
    0.42 0.67 0.60
    0.30 0.49 0.74
    0.47 0.63 0.67
    0.36 0.55 0.44
];
for k = 1:m
    plot(0:T-1, meanA(:,k), 'LineWidth', 2, 'Color', colorsA(k,:));
end
for k = 1:m
    plot(0:T-1, meanD(:,k), 'LineWidth', 2, 'Color', colorsD(k,:));
end
legend_entries = cell(1,2*m);
for k = 1:m
    legend_entries{k} = ['Adopters [', num2str(k), ']'];
end
for k = 1:m
    legend_entries{m+k} = ['Dissatisfied [', num2str(k), ']'];
end
legend(legend_entries, 'Location', 'best', 'NumColumns',2);
xlabel('Time', 'Interpreter','latex'); ylabel('Fraction of population','Interpreter','latex'); 
grid on; set(gca,'FontSize',15); ylim([0 0.7]);
