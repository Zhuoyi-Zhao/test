function optimize_weighted_AoI
    % Given parameters (example):
    N = 5; 
    p = [0.8; 0.4; 0.9; 0.7; 0.6];    % ensure column vector
    L = [5; 3; 4; 6; 2];              % ensure column vector
    w = [1; 2; 1; 1.5; 0.5];          % ensure column vector
    
    % Compute E[Y_j^2] for all j:
    E_Y2 = zeros(N,1);
    for j = 1:N
        E_Y2(j) = 1 + 2*(L(j)-1) + ((L(j)-1)*(L(j)-p(j)))/p(j);
    end
    
    % Optimization variables: mu
    lb = zeros(N,1);
    ub = ones(N,1);
    Aeq = ones(1,N);
    beq = 1;
    mu0 = ones(N,1)/N; % initial guess
    
    % Objective function
    objfun = @(mu) weighted_AoI(mu, N, p, L, E_Y2, w);
    
    options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
    [mu_opt, fval_opt] = fmincon(objfun, mu0, [], [], Aeq, beq, lb, ub, [], options);
    
    fprintf('Optimal mu:\n');
    disp(mu_opt');
    fprintf('Optimal weighted average AoI: %f\n', fval_opt);
    
    % Compute AoI for each source at the optimal mu
    E_h_opt = zeros(1,N);
    for i = 1:N
        E_h_opt(i) = E_h_i(mu_opt, i, p, L, E_Y2);
    end
    weighted_avg_AoI = (1/N)*sum(w.*E_h_opt');
    fprintf('Check Weighted Avg AoI (direct computation): %f\n', weighted_avg_AoI);
end

function val = weighted_AoI(mu, N, p, L, E_Y2, w)
    E_h = zeros(N,1);
    for i = 1:N
        E_h(i) = E_h_i(mu, i, p, L, E_Y2);
    end
    val = (1/N)*sum(w.*E_h);
end

function val = E_h_i(mu, i, p, L, E_Y2)
    % Compute E[h_i(t)] using the derived formula.

    % U = sum_j mu_j L_j
    U = sum(mu.*L);
    
    % A_i and B_i:
    A_i = (U - mu(i)*(L(i)-1)) / (mu(i)*p(i));
    B_i = (L(i)-1)/p(i);
    
    E_Si2 = ((L(i)-1)*(L(i)-p(i)))/p(i)^2;
    
    % Compute E[(Wi[m])^2]:
    % sum_{j != i} mu_j(E[Y_j^2]+2L_j A_i)
    N = length(p);
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