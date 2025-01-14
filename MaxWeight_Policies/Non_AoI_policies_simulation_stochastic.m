function [actual_throughput,optimal_cost]=Non_AoI_policies_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_simplified,set_Policy)
%[actual_throughput,optimal_cost]=MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_simplified,set_Policy)
%set_Policy=1 for ZERO DELAY policy with maximum throughput, set_Policy=0 for Minimum Delay policy
%set_simplified chooses between updating the AoI at every slot or ONLY at the beginning of frames.

%Author: Igor Kadota
%Advisor: Eytan Modiano
%Period: Fall 2018

global probability cost alpha arrival

%% Simulation Setup
%num_iterations=100;
%set_simplified=0; %set_simplified=1 --> AoI is updated at the end of frames
%set_Policy=1 --> Max Weight
%set_Whittle=0 --> Greedy on p*Debt


%% Model Setup
%K=1000;
%T=10;
%h1=1;
%M=2;
%age_weight=1;
%debt_weight=1;

p=zeros(M,1);
p=[p+probability;0]; %Making sure that "probability" has the right size

A=zeros(M,1);
A =[A+alpha;0]; %Making sure that "alpha" has the right size

C=zeros(M,1);
C =[C+cost;0]; %Making sure that "cost" has the right size

Arr=zeros(M,1);
Arr =[Arr+arrival;0]; %Making sure that "arrival" has the right size

optimal_cost=0;

%% Simulation

real_q=zeros(M,num_iterations);
AoI=zeros(M,num_iterations);
node_index=0;

for iteration=1:1:num_iterations

    %transmissions=zeros(M+1,1);
    delivered=zeros(M+1,1);    
    throughput_debt=zeros(M,1);
    h=h1*ones(M,1);
    z=z1*ones(M,1);
    vector_successful_tx=zeros(M,1);

    for frame_index=1:1:K       
        %vector_delayed_ACK=(M+1)*ones(T,1);
        
        for time=1:1:T            
            for node_index=1:1:M
                if rand<Arr(node_index)
                    z(node_index)=0;
                    vector_successful_tx(node_index)=0;% new arrival
                else
                    z(node_index)=z(node_index)+1;
                end
            end
            if sum(vector_successful_tx)==M
            %Idle
                node_index=M+1;                                    
            elseif set_Policy==1 %The ZERO DELAY policy with maximum throughput
                rank=(A(1:M).*p(1:M)).*(z==0).*(1-vector_successful_tx);                                
                if max(rank)==0
                    node_index=M+1;
                else
                    possible_nodes=find(rank==max(rank));
                    node_index=possible_nodes(ceil(rand*length(possible_nodes)));
                end
                
                %[tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*(h-z).*(1-vector_successful_tx));
            elseif set_Policy==0 %The Minimum Delay policy
                %if sum(max(throughput_debt(1:M),0))>0
                rank=z.*(1-vector_successful_tx);
                possible_nodes=find(rank==min(rank));
                node_index=possible_nodes(ceil(rand*length(possible_nodes)));
                %else
                %    node_index=mod(node_index,M)+1; %Do something similar to Round Robin
                %end            
            else
                disp('There is a problem in set_Policy. set_Whittle=1 for Max Weight, set_Whittle=0 for Greedy on P*Debt.')
                break;
            end
            
            if node_index~=(M+1)
                AoI(node_index,iteration)=AoI(node_index,iteration)+C(node_index);
            end
            
            if set_simplified==0
            %AoI is updated at every SLOT
            if node_index~=(M+1)
                AoI((vector_successful_tx(:,1)==0),iteration)=AoI((vector_successful_tx(:,1)==0),iteration)+A((vector_successful_tx(:,1)==0)).*h((vector_successful_tx(:,1)==0));
            end
            else
            %AoI is updated at the end of every FRAME
                AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h(:);
            end          
            
            %transmissions(node_index_greedy,1)=transmissions(node_index_greedy,1)+1;            
            
            if (node_index~=(M+1)) && (rand<p(node_index)) && (vector_successful_tx(node_index)==0)
                delivered(node_index)=delivered(node_index)+1;
                vector_successful_tx(node_index)=1;
            end  
        end
        h(vector_successful_tx(:,1)==1)=z(vector_successful_tx(:,1)==1)+1;
        h(vector_successful_tx(:,1)==0)=h(vector_successful_tx(:,1)==0)+1;                
    end
    real_q(:,iteration)=delivered(1:M)/K;
end
AoI_average_GRE=mean(AoI,2); %Averaging on the ITERATIONS
optimal_cost_clients=T*AoI_average_GRE(:)/K;
optimal_cost=sum(optimal_cost_clients);
actual_throughput=mean(real_q,2);
%optimal_cost_GRE=sum(q_vector_GRE(:,1).*alpha(1:M));

end