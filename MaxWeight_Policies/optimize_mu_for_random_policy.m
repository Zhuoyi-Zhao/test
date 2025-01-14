function mu_opt = optimize_mu_for_random_policy(N, p, L, alpha, E_Y2)
% This function optimizes mu to minimize the weighted average AoI 
% under the no-switch random policy scenario.
% Inputs:
%   N, p, L, alpha, E_Y2
% Output:
%   mu_opt: optimal probability distribution (Nx1)

lb = zeros(N,1);
ub = ones(N,1);
Aeq = ones(1,N);
beq = 1;
mu0 = ones(N,1)/N;

objfun = @(mu) weighted_AoI(mu, p, L, E_Y2, alpha);

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
mu_opt = fmincon(objfun, mu0, [], [], Aeq, beq, lb, ub, [], options);
end



