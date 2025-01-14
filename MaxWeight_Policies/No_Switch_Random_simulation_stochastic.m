function [actual_throughput,optimal_cost]=No_Switch_Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,mu_opt_no_switch)
% This function simulates the NO-SWITCH RANDOMIZED policy.
%
% Inputs and outputs are the same as the original function provided by user.
% The difference is in the policy:
% - We only choose a node when no packet is in service.
% - Once chosen, if the first transmission fails, next slot we choose again.
% - If first transmission succeeds, we continue with the same node until
%   all L_i packets are completed (no switching in between).
%
% Global variables:
% probability: p_i for each source i
% alpha: weight (already given)
% arrival: arrival rates
% optMU: the mu vector (probabilities to choose a node)
% L: packets per arrival

global probability cost alpha arrival optQ  L optMU

% Setup
p=zeros(M,1);
p=[p+probability;0]; 
A=zeros(M,1);
A=[A+alpha;0]; 
C=zeros(M,1);
C=[C+cost;0]; 
Arr=zeros(M,1);
Arr=[Arr+arrival;0];

scheduling_probabilities=mu_opt_no_switch;

optimal_cost=0;

real_q=zeros(M,num_iterations);
AoI=zeros(M,num_iterations);

for iteration=1:num_iterations
    delivered=zeros(M+1,1);    
    h=h1*ones(M,1); % AoI per source
    z=zeros(M,1);
    vector_successful_tx=zeros(M,1); % 0 means no ongoing packet
    current_node = M+1; % M+1 means idle, no packet chosen yet

    for slot_index=1:K
        % Arrival: if no ongoing packet for a node and arrival occurs, start a new L(i) packet:


        % Choosing node logic for NO-SWITCH:
        % If current_node = M+1 or current_node finished last slot (no ongoing packet),
        % we pick a new node based on mu if available
        if current_node==(M+1) || vector_successful_tx(current_node)==0
                % pick a node according to mu
            randomizer=rand*sum(scheduling_probabilities);
            idx=0; rand_value=0;
            while randomizer>rand_value
                idx=idx+1;
                rand_value=rand_value+scheduling_probabilities(idx);
            end
            current_node=idx;
        end

        for node_index_Renew=1:M
            if rand<Arr(node_index_Renew) && vector_successful_tx(node_index_Renew)==0
                vector_successful_tx(node_index_Renew)=L(node_index_Renew);
            end
        end

        % AoI update:
        AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h(:);
        h(:)=h(:)+1; % age increments by 1 each slot

        % Transmission attempt if current_node != M+1
        if (current_node~=(M+1)) && (vector_successful_tx(current_node)==L(current_node))
            % Attempt a transmission:
            if rand<p(current_node)
                % success
                vector_successful_tx(current_node)=vector_successful_tx(current_node)-1;
                if vector_successful_tx(current_node)==0
                    % Packet completed:
                    h(current_node)=z(current_node)+1;
                    z(current_node)=0;
                    delivered(current_node)=delivered(current_node)+1;
                    % Next slot we can choose again
                else
                    % Still have packets for the same node, no switching,
                    % continue to next slot with same node
                end
            else
                current_node=M+1;
            end
        elseif (current_node~=(M+1)) && (vector_successful_tx(current_node)>0) && (vector_successful_tx(current_node)<L(current_node))
            if rand<p(current_node)
                % success
                vector_successful_tx(current_node)=vector_successful_tx(current_node)-1;
                if vector_successful_tx(current_node)==0
                    % Packet completed:
                    h(current_node)=z(current_node)+1;
                    z(current_node)=0;
                    delivered(current_node)=delivered(current_node)+1;
                    % Next slot we can choose again
                else
                    % Still have packets for the same node, no switching,
                    % continue to next slot with same node
                end
            end
        end

        % Update z for partially completed:
        for nd=1:M
            if vector_successful_tx(nd)<L(nd) && vector_successful_tx(nd)>0
                z(nd)=z(nd)+1;
            end
        end
    end
    real_q(:,iteration)=delivered(1:M)/K;
end

AoI_average_GRE=mean(AoI,2);
optimal_cost_clients=T*AoI_average_GRE(:)/K;
optimal_cost=sum(optimal_cost_clients)/M;
actual_throughput=mean(real_q,2);

end