function [actual_throughput,optimal_cost]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations)
%[actual_throughput,optimal_cost]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,work_conserving_policy)
%work_conserving_policy=1 then the randomized policy is work conserving and only attempts transmission of UNDELIVERED packets

%Author: Igor Kadota
%Advisor: Eytan Modiano
%Period: Fall 2018

global probability cost alpha arrival optMU L

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

scheduling_probabilities=optMU;

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
    z=zeros(M,1);
    vector_successful_tx=L;

    for slot_index=1:1:K       
                
        for node_index=1:1:M
            
            if rand<Arr(node_index) && vector_successful_tx(node_index)==0 
                % z(node_index)=0;
                vector_successful_tx(node_index)=L(node_index);% new arrival
            % else 
                % z(node_index)=z(node_index)+1;
            end
        end
        if sum(vector_successful_tx)==0 %Idle
            node_index=M+1;                                    
        else %The Randomized Policy is used                
            randomizer=rand*sum(scheduling_probabilities);
            node_index=0;
            rand_value=0;
            while randomizer>rand_value
                node_index=node_index+1;
                rand_value=rand_value+scheduling_probabilities(node_index);
            end             
            %if the client was already transmitted in this frame, I choose another client.
            % while (vector_successful_tx(node_index)==0)&&(work_conserving_policy==1)
            %     randomizer=rand*sum(scheduling_probabilities);
            %     node_index=0;
            %     rand_value=0;
            %     while randomizer>rand_value
            %         node_index=node_index+1;
            %         rand_value=rand_value+scheduling_probabilities(node_index);
            %     end
            % end
        end       
    
        AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h(:);       
        h(:)=h(:)+1;                
        
        if (node_index~=(M+1)) && (rand<p(node_index)) %&& (vector_successful_tx(node_index)>=1)
            vector_successful_tx(node_index)=vector_successful_tx(node_index)-1;
            if vector_successful_tx(node_index)==0
                h(node_index)=z(node_index)+1;
                z(node_index)=0;
                delivered(node_index)=delivered(node_index)+1;
            end
        end  
        for node_index = 1:1:M
            if vector_successful_tx(node_index)<L(node_index) && vector_successful_tx(node_index)>0
                z(node_index) = z(node_index)+1;
            end
        end
    end
    real_q(:,iteration)=delivered(1:M)/K;
end
AoI_average_GRE=mean(AoI,2); %Averaging on the ITERATIONS
optimal_cost_clients=T*AoI_average_GRE(:)/K;
optimal_cost=sum(optimal_cost_clients);
% optimal_cost_Rand=(sum(sqrt(A(1:M)./p(1:M).*(2*L.^2+L-1)./(2.*L))))^2+sum(A(1:M).*(1./Arr(1:M)-1));
actual_throughput=mean(real_q,2);
%optimal_cost_GRE=sum(q_vector_GRE(:,1).*alpha(1:M));

end