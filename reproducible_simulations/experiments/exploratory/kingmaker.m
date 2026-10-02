%%%%%%%%%%%%%%%%%%%% ADOZIONE m>=2 tecnologie
%%%%%%%% KING MAKER EFFECT
%% Parametri
rng(3);
n = 30;          % numero di agenti
m = 3;           % numero di tecnologie
T = 600;         % orizzonte temporale

% entry time di ciascuna tecnologia
% es: tech 1 presente da subito, tech 2 entra a t=100, tech 3 a t=150, tech 4 a t=200
% t_entry = [1, 100, 150, 200];

%% Parametri per ogni tecnologia (n x m)
beta   = rand(n,m); 
beta(:,1) = 0.3+0.2*rand(n,1);
beta(:,2) = 0.4+0.2*rand(n,1);
delta  = rand(n,m);
delta(:,1) = 0.1+ 0.3*rand(n,1);
delta(:,2) = 0.2 + 0.4*rand(n,1);
gamma  = 0.2*rand(n,m);
gamma(:,2) = 0.3*rand(n,1);
gamma(:,3) = 0.6*rand(n,1);

lambda = rand(n,m);
lambda = 0.8*lambda ./ sum(lambda,2);   % ogni riga somma a 1
xi = rand(n,m) .* (1 - lambda);
% 
%% Matrici di rete
W = rand(n,n);
W = W ./ sum(W,2);

tilde_W = rand(n,n);
tilde_W = tilde_W ./ sum(tilde_W,2);

% Stati
% S(t,i)         = suscettibili/non adottanti
% A(t,i,k)       = adottanti tecnologia k
% D(t,i,k)       = insoddisfatti della tecnologia k
% X(t,i,k)       = opinione sulla tecnologia k

S = zeros(T,n);
A = zeros(T,n,m);
D = zeros(T,n,m);
X = zeros(T,n,m);

%% Inizializzazione
% for k = 1:m
%     A(1,:,k) = 0.15*rand(1,n);
% end
A(1,:,1) = rand(1,n);
A(1,:,2) = 0.5*rand(1,n);
A(1,:,3) = zeros(1,n);
% A(1,:,4) = zeros(1,n);
% A(1,:,5) = 0.1*rand(1,n);
S(1,:) = max(zeros(1,n),ones(1,n)-sum(A(1,:,:),3));
D(1,:,:) = 0;
X(1,:,:) = 0.4;

%% Simulazione
for t = 1:T-1
    
    % --- calcolo forze di adozione per tutte le tecnologie ---
    % Wa(i,k) = sum_j W(i,j) * A(t,j,k)
    Wa = zeros(n,m);
    for k = 1:m
        Wa(:,k) = W * squeeze(A(t,:,k))';
    end
    
    for i = 1:n
        s_old = S(t,i);
        adoption_force = 0;
        if t <= 300
            for k = 1:2
                %% 1) aggiornamento opinioni
                X0ik = X(1,i,k);
                X_neighbors = tilde_W(i,:) * squeeze(X(t,:,k))';
                X(t+1,i,k) = (1 - lambda(i,k) - xi(i,k)) * X0ik + lambda(i,k) * X_neighbors + xi(i,k) * Wa(i,k);
                X(t+1,i,k) = max(0, min(1, X(t+1,i,k)));
    
                %% 2) aggiornamento non-adottanti
                adoption_force = adoption_force + beta(i,k) * X(t,i,k) * Wa(i,k);
                S(t+1,i) = s_old - s_old * adoption_force;
                S(t+1,i) = max(0, min(1, S(t+1,i)));
            
                %% 3) aggiornamento adottanti A^{[k]}
                dissatisfied_others = 0;
                for h = 1:m
                    if h ~= k
                        dissatisfied_others = dissatisfied_others + D(t,i,h);
                    end
                end
    
                A(t+1,i,k) = A(t,i,k) + beta(i,k) * X(t,i,k) * s_old * Wa(i,k) - delta(i,k) * A(t,i,k) + gamma(i,k) * X(t,i,k) * dissatisfied_others;
                A(t+1,i,k) = max(0, min(1, A(t+1,i,k)));
    
                %% 4) aggiornamento insoddisfatti D^{[k]}
                % tasso con cui un insoddisfatto di k passa ad altre tecnologie
                switching_out = 0;
                for h = 1:m
                    if h ~= k
                        switching_out = switching_out + gamma(i,h) * X(t,i,h);
                    end
                end
              
                D(t+1,i,k) = D(t,i,k) + delta(i,k) * A(t,i,k) - D(t,i,k) * switching_out;
                D(t+1,i,k) = max(0, min(1, D(t+1,i,k)));
            end
        else
            for k = 1:m
                %% 1) aggiornamento opinioni
                X0ik = X(1,i,k);
                X_neighbors = tilde_W(i,:) * squeeze(X(t,:,k))';
                X(t+1,i,k) = (1 - lambda(i,k) - xi(i,k)) * X0ik + lambda(i,k) * X_neighbors + xi(i,k) * Wa(i,k);
                X(t+1,i,k) = max(0, min(1, X(t+1,i,k)));
    
                %% 2) aggiornamento non-adottanti
                adoption_force = adoption_force + beta(i,k) * X(t,i,k) * Wa(i,k);
                S(t+1,i) = s_old - s_old * adoption_force;
                S(t+1,i) = max(0, min(1, S(t+1,i)));
            
                %% 3) aggiornamento adottanti A^{[k]}
                dissatisfied_others = 0;
                for h = 1:m
                    if h ~= k
                        dissatisfied_others = dissatisfied_others + D(t,i,h);
                    end
                end
    
                A(t+1,i,k) = A(t,i,k) + beta(i,k) * X(t,i,k) * s_old * Wa(i,k) - delta(i,k) * A(t,i,k) + gamma(i,k) * X(t,i,k) * dissatisfied_others;
                A(t+1,i,k) = max(0, min(1, A(t+1,i,k)));
    
                %% 4) aggiornamento insoddisfatti D^{[k]}
                % tasso con cui un insoddisfatto di k passa ad altre tecnologie
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
end

%% Medie
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

% Plot adottanti + insoddisfatti
figure; hold on;
% Palette calda per Adopters
colorsA = [
    0.75 0.33 0.25   % terracotta
    0.86 0.49 0.20   % burnt orange
    0.93 0.69 0.13   % warm gold
    0.80 0.52 0.25   % amber brown
    0.91 0.57 0.33   % soft orange
    0.70 0.42 0.18   % copper
];
% Palette fredda per Dissatisfied
colorsD = [
    0.16 0.44 0.52   % deep teal
    0.24 0.60 0.56   % jade teal
    0.42 0.67 0.60   % soft green
    0.30 0.49 0.74   % muted blue
    0.47 0.63 0.67   % dusty cyan
    0.36 0.55 0.44   % sage green
];
for k = 1:m
    plot(0:T-1, meanA(:,k), 'LineWidth', 2, 'Color', colorsA(k,:));
end
for k = 1:m
    plot(0:T-1, meanD(:,k), 'LineWidth', 2, 'Color', colorsD(k,:));
end
% Plot susceptibles
%plot(0:T-1, meanS, 'k', 'LineWidth', 2);
legend_entries = cell(1,2*m);
for k = 1:m
    legend_entries{k} = ['Adopters [', num2str(k), ']'];
end
for k = 1:m
    legend_entries{m+k} = ['Dissatisfied [', num2str(k), ']'];
end
%legend_entries{2*m+1} = 'Susceptibles';
legend(legend_entries, 'Location', 'best', 'NumColumns',2);
xlabel('Time', 'Interpreter','latex'); ylabel('Fraction of population','Interpreter','latex'); 
grid on; set(gca,'FontSize',15); ylim([0 0.7]);
%saveas(gcf, 'adopters4.png')

% %% Plot opinioni
% figure; hold on;
% for k = 1:m
%     plot(0:T-1, meanX(:,k), 'LineWidth', 2);
% end
% xlabel('Time','Interpreter','latex');
% %ylabel('Average opinion');
% legend(arrayfun(@(k) ['Opinion [', num2str(k), ']'], 1:m, 'UniformOutput', false), ...
%        'Location', 'best');
% grid on; set(gca,'FontSize',15); %saveas(gcf, 'opinions4.png')