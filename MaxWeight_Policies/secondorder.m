clear; clc; rng('shuffle');

%% Parameters
N = 3;               % number of sources
j = 1;               % target source index
mu = [0.4, 0.3, 0.3];% selection probabilities
p = [0.8, 0.6, 0.9]; % success probabilities
L = [5, 3, 8];       % number of packets for each source

% Verify that sum(mu)=1
assert(abs(sum(mu)-1)<1e-12,'Mu must sum to 1');

numTrials = 6e6;    % number of Monte Carlo trials for simulation

%% Theoretical Computation (Revised Formula)

% First moment:
% E[S_j] = (mu_j + sum_{i!=j} mu_i L_i)/(mu_j p_j)
sum_L = sum(mu((1:N)~=j).*L((1:N)~=j));
E_Sj = (mu(j) + sum_L)/(mu(j)*p(j));

% Compute E[Y_i^2] for each i
EYi2 = zeros(1,N);
for i = 1:N
    if i==j
        % For the target j, Y_j is not defined the same way, we don't need E[Y_j^2].
        EYi2(i) = NaN; 
    else
        % E[Y_i^2] = 1 + 2(L_i-1) + ((L_i-1)(L_i - p_i))/p_i
        EYi2(i) = 1 + 2*(L(i)-1) + ((L(i)-1)*(L(i)-p(i)))/p(i);
    end
end

% Compute the second moment using revised formula:
% E[S_j^2] = 2E[S_j]^2 - 2E[S_j] + (mu_j + sum_{i!=j} mu_i E[Y_i^2])/(mu_j p_j)

sum_EYi2 = sum(mu((1:N)~=j).*EYi2((1:N)~=j));
E_Sj2 = 2*E_Sj^2 - 2*E_Sj + (mu(j) + sum_EYi2)/(mu(j)*p(j));

%% Simulation
S_values = zeros(numTrials,1);

for trial = 1:numTrials
    S = 0;
    done = false;
    while ~done
        % Choose a source
        i = find(rand<=cumsum(mu),1,'first');
        
        if i == j
            % Attempt once
            S = S + 1;
            if rand < p(j)
                % success
                done = true;
            else
                % fail, continue
            end
            
        else
            % Non-target source i
            % First attempt:
            S = S + 1;
            if rand < p(i)
                % Success first packet, now transmit L_i-1 packets
                remaining = L(i)-1;
                for rr = 1:remaining
                    % Each packet: geometric trials until success
                    while true
                        S = S+1;
                        if rand < p(i)
                            break; 
                        end
                    end
                end
                % After finishing source i, back to selection
            else
                % failed the first attempt of source i
                % back to selection immediately
            end
        end
    end
    S_values(trial) = S;
end

sim_mean = mean(S_values);
sim_mean2 = mean(S_values.^2);

%% Results Comparison
fprintf('Theoretical E[S_j]: %f\n', E_Sj);
fprintf('Simulated  E[S_j]: %f\n', sim_mean);
fprintf('Difference: %f\n', abs(E_Sj - sim_mean));

fprintf('\nTheoretical E[S_j^2]: %f\n', E_Sj2);
fprintf('Simulated  E[S_j^2]: %f\n', sim_mean2);
fprintf('Difference: %f\n', abs(E_Sj2 - sim_mean2));