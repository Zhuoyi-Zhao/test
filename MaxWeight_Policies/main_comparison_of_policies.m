function main_comparison_of_policies
    % This main function compares various policies (maxweight, DPP, Age-debt, etc.)
    % with the no-switch randomized policy we have derived. We will:
    % 1. Define system parameters.
    % 2. Compute the E[Y_j^2].
    % 3. Compute or load results from other policies as in previous code.
    % 4. Call our optimize_mu_for_random_policy to get the optimal mu for the no-switch random strategy.
    % 5. Compute weighted average AoI for that mu and add it to the comparison plot.

    clear; clc;
    global L
    
    % Example system parameters:
    N = 4;
    p = [0.8; 0.4; 0.9; 0.7];     % success probabilities (as column vector)
    L = [3; 5; 7; 4];             % packet counts
    w = [4;8;2;1];          % weights for AoI
    
    % Compute E[Y_j^2] for all j
    E_Y2 = zeros(N,1);
    for j=1:N
        E_Y2(j) = 1 + 2*(L(j)-1) + ((L(j)-1)*(L(j)-p(j)))/p(j);
    end
    
    % Suppose we have results for other policies already computed (just placeholders):
    vector_p = 0.05:0.05:0.5; 
    optimal_cost_OtherPolicy = rand(1,length(vector_p))*5 + 10; % placeholder for some other policy results

    % Now integrate the random no-switch policy optimization:
    % For each p_1 in vector_p, we do:
    optimal_cost_RandomNoSwitch = zeros(1,length(vector_p));
    for idx = 1:length(vector_p)
        p(1) = vector_p(idx); % vary p(1)
        
        % Optimize mu using previously derived approach:
        mu_opt = optimize_mu_for_random_policy(p, L, w, E_Y2);
        
        % Once we have mu_opt, compute weighted average AoI:
        AoI = weighted_AoI(mu_opt, p, L, E_Y2, w);
        optimal_cost_RandomNoSwitch(idx) = AoI;
    end
    
    % Plot comparison:
    figure;
    hold on;
    plot(vector_p, optimal_cost_OtherPolicy, 'r-o', 'LineWidth',2,'MarkerSize',10);
    plot(vector_p, optimal_cost_RandomNoSwitch, 'b--^', 'LineWidth',2,'MarkerSize',10);
    legend('Other Policy','Random No-Switch Optimal','Location','Best');
    xlabel('p_1');
    ylabel('Weighted Avg AoI');
    title('Comparison of policies');
    hold off;

end

function mu_opt = optimize_mu_for_random_policy(p, L, w, E_Y2)
    % This function uses fmincon to find the optimal mu that minimizes 
    % the weighted average AoI for the no-switch random policy scenario.
    %
    % Inputs:
    %   p, L, w, E_Y2: system parameters
    % Output:
    %   mu_opt: optimal probability vector of choosing each source
    
    N = length(p);
    lb = zeros(N,1);
    ub = ones(N,1);
    Aeq = ones(1,N);
    beq = 1;
    mu0 = ones(N,1)/N; % initial guess
    
    % Objective function: Weighted AoI given mu
    objfun = @(mu) weighted_AoI(mu, p, L, E_Y2, w);
    
    options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
    mu_opt = fmincon(objfun, mu0, [], [], Aeq, beq, lb, ub, [], options);
end

function val = weighted_AoI(mu, p, L, E_Y2, w)
    % Compute (1/N)*sum w_i E[h_i(t)] given mu
    N = length(p);
    E_h = zeros(N,1);
    for i = 1:N
        E_h(i) = E_h_i(mu, i, p, L, E_Y2);
    end
    val = (1/N)*sum(w.*E_h);
end

function val = E_h_i(mu, i, p, L, E_Y2)
    % Compute E[h_i(t)] using previously derived formula
    % We have (as from previous derivations):
    % E[h_i(t)] = (mu_i p_i / U)*[ ((E[Wi^2]+E[Si^2])/2) + B_i^2 + 2 A_i B_i ] + 1
    %
    % E[Si[m]] = (L_i-1)/p_i = B_i
    % E[(Si[m])^2] = ((L_i -1)(L_i-p_i))/p_i^2
    %
    % E[(Wi[m])^2] = [ mu_i(1 + 2(1-p_i)A_i) + sum_{j!=i} mu_j(E[Y_j^2]+2L_j A_i) ] / (mu_i p_i)
    %
    % A_i = (U - mu_i(L_i-1)) / (mu_i p_i)
    % B_i = (L_i -1)/p_i
    % U = sum_j mu_j L_j

    N = length(p);
    U = sum(mu.*L);
    A_i = (U - mu(i)*(L(i)-1)) / (mu(i)*p(i));
    B_i = (L(i)-1)/p(i);
    
    E_Si2 = ((L(i)-1)*(L(i)-p(i)))/p(i)^2;
    
    % sum_part for E[Wi^2]:
    sum_part = 0;
    for jj = 1:N
        if jj ~= i
            sum_part = sum_part + mu(jj)*(E_Y2(jj) + 2*L(jj)*A_i);
        end
    end
    E_Wi2 = ( mu(i)*(1 + 2*(1-p(i))*A_i) + sum_part )/(mu(i)*p(i));
    
    numerator_part = ((E_Wi2 + E_Si2)/2) + B_i^2 + 2*A_i*B_i;
    E_hi = (mu(i)*p(i)/U)*numerator_part + 1;
    val = E_hi;
end