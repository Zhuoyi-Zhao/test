clear
clc

global probability cost alpha arrival optQ optMU optMU_buff optMU_NO_buff L

tic;

%% Model Setup
M = 4;
vector_p = 0.05:0.05:1;
L = [3;5;7;4];
num_setups = length(vector_p);
num_iterations = 10;
T = 1;
set_simplified = 1;
h1 = 1; % Initial age
z1 = 0; % MUST BE LOWER THAN H1

%% Parameters Setup
optimal_cost_Rand = zeros(1, num_setups);
optimal_cost_Max_Weight_Age_debt = zeros(1, num_setups);
optimal_cost_LowerBound_stoch = zeros(1, num_setups);

%% Simulations
for count = 1:num_setups
    count
    K = 100000;
    
    cost = zeros(M,1);
    alpha=[4;8;2;1];
    arrival=[1;1;1;1];
    probability=[vector_p(count);0.4;0.6;0.8]; 
    p = probability;
    A = alpha;
    Arr = arrival;
    
    Lower_Bound_find_q(M); % Finding the optQ
    optimal_cost_LowerBound_stoch(count) = sum(alpha .* (1 + 1 ./ optQ)) / (2 * M);
    disp('Lower Bound for Stochastic Arrivals')
    
    optMU = sqrt(A ./ p .* (3 * L.^2 - L) ./ (2 .* L));
    sum_optMU = sum(optMU);
    optMU = optMU / sum_optMU;
    
    optMU_NO_buff = sqrt(alpha ./ (arrival .* probability));
    sum_optMU_NO_buff = sum(optMU_NO_buff);
    optMU_NO_buff = optMU_NO_buff / sum_optMU_NO_buff;
    
    Random_find_mu_buffered_BerBer1(M) % Finding the optMU_buff
    optimal_cost_Rand(count) = (sum(sqrt(A ./ p .* (3 * L.^2 - L) ./ (2 .* L))))^2 / M + sum(A .* 3 / 2) / M;
    disp('Randomized Policy')
    
    %% Optimize Achievable_AoI using Gradient Descent
    set_Policy = 0; % 0: Age-Debt Policy
    % Initial Achievable_AoI
    Achievable_AoI = (1 ./ p ./ optMU .* (3 * L.^2 - L) ./ (2 .* L))/5 + 1 * 3 / 2;
    
    % Set gradient descent parameters
    max_iter = 1000;
    tolerance = 1e-6;
    learning_step = Achievable_AoI./Achievable_AoI; % Learning rate, may need to adjust for convergence
    previous_cost = 0;
    
    for iter = 1:max_iter
        % Run simulation to get current per_stream_AoI and total_cost
        set_Policy = 0;
        [Q, total_cost] = LQ_MaxWeight_and_Greedy_simulation_stochastic(K, T, M, h1, z1, num_iterations, set_Policy, Achievable_AoI);
        total_cost = total_cost / (M * T);
        
        % Check if total_cost is NaN or Inf
        if isnan(total_cost) || isinf(total_cost)
            warning('Total cost is NaN or Inf, adjusting learning rate or initial Achievable_AoI may help.');
            break;
        end
        
        % Calculate the gradient (derivative of the loss function with respect to Achievable_AoI)
        % Assuming the loss function is total_cost = sum(alpha .* per_stream_AoI)
        % The gradient is simply alpha
        % gradient = (per_stream_AoI - Achievable_AoI);
        % gradient = Achievable_AoI - per_stream_AoI;                                             
        
        % Update Achievable_AoI
        for node_index = 1:1:M   
            if Q(node_index)>Achievable_AoI(node_index)*20
                Achievable_AoI(node_index) = Achievable_AoI(node_index) +  learning_step(node_index);
            else
                Achievable_AoI(node_index) = Achievable_AoI(node_index) -  learning_step(node_index);
                learning_step = learning_step * 0.9;
                disp('Total cost increased, reducing learning rate.'); 
            end
        end
        fprintf('Iteration %d, Total Cost: %f\n ', iter, total_cost);
        % Ensure Achievable_AoI is non-negative
        Achievable_AoI = max(Achievable_AoI, 0)
        
        % Print information about the current iteration

        
        % % Check for convergence
        % if abs(previous_cost - total_cost) < tolerance
        %     disp('Converged.');
        %     break;
        % end
        
        % If total_cost increased, reduce the learning rate
        % if total_cost > previous_cost
        %     % learning_rate = learning_rate * 0.9;
        %     disp('Total cost increased, reducing learning rate.');
        % end
        
        previous_cost = total_cost;
    end
    
    % Store the final cost
    optimal_cost_Max_Weight_Age_debt(count) = total_cost;
    disp(['Age-Debt Policy optimized with Gradient Descent, final cost: ', num2str(total_cost)]);
    toc
end

% Plotting
figure(4)
hold on
[Rand_plot] = plot(vector_p, optimal_cost_Rand, 'bo--', 'LineWidth', 3, 'MarkerSize', 10);
[Age_Debt_plot] = plot(vector_p, optimal_cost_Max_Weight_Age_debt, 'rx--', 'LineWidth', 3, 'MarkerSize', 10);
[Lower_Bound] = plot(vector_p, optimal_cost_LowerBound_stoch, 'k-', 'LineWidth', 3);
legend([Rand_plot, Age_Debt_plot, Lower_Bound], 'Randomized', 'Age-Debt (Optimized)', 'Lower Bound', 'Location', 'NorthWest');
ylabel('Expected Weighted Sum AoI')
xlabel('Channel Reliability of Stream 1, p_1')
hold off

toc;