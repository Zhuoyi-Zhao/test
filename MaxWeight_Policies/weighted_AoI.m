function val = weighted_AoI(mu, p, L, E_Y2, alpha)
% Compute (1/N)*sum alpha_i * E[h_i(t)] given mu
N = length(p);
E_h = zeros(N,1);
for i = 1:N
    E_h(i) = E_h_i(mu, i, p, L, E_Y2);
end
val = (1/N)*sum(alpha.*E_h);
end
