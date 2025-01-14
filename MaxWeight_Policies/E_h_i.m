function val = E_h_i(mu, i, p, L, E_Y2)
% Compute E[h_i(t)] using previously derived formula:
% E[h_i(t)] = (mu_i p_i / U)*[ ((E[Wi^2]+E[Si^2])/2) + B_i^2 + 2 A_i B_i ] + 1
N = length(p);
U = sum(mu.*L);
A_i = (U - mu(i)*(L(i)-1)) / (mu(i)*p(i));
B_i = (L(i)-1)/p(i);

E_Si2 = ((L(i)-1)*(L(i)-p(i)))/p(i)^2;

% sum_part for E[Wi^2]
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