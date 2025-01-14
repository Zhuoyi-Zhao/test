function [actual_throughput,optimal_cost]=MaxWeight_Buffered_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy)
%[actual_throughput,optimal_cost]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,work_conserving_policy)
%set_Policy=1 for LINEAR Max_Weight of AoI with buffers, set_Policy=0 for Throughput MaxWeight (send the one with longest queue)

%Author: Igor Kadota
%Advisor: Eytan Modiano
%Period: Fall 2018

global probability alpha arrival optMU_buff

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

Arr=zeros(M,1);
Arr =[Arr+arrival;0]; %Making sure that "arrival" has the right size

optimal_cost=0;

%% Simulation

real_q=zeros(M,num_iterations);
AoI=zeros(M,num_iterations);
node_index=0;

for iteration=1:1:num_iterations

    %transmissions=zeros(M+1,1);
    delivered=zeros(M,1);
    throughput_debt=zeros(M,1);
    h=h1*ones(M,1);
    z=z1*ones(M,1); %Delay of the FRESHEST packet at the source (a.k.a. AoI of the source)
    %vector_successful_tx=zeros(M,1);
    arrived=zeros(M,1);
    delay_buffer=zeros(M,K);

    for slot_index=1:1:K       
                
        for node_index=1:1:M
            if rand<Arr(node_index)
                z(node_index)=0;
                %vector_successful_tx(node_index)=0;% new arrival
                arrived(node_index)=arrived(node_index)+1;
                delay_buffer(node_index,arrived(node_index))=slot_index;
            else
                z(node_index)=z(node_index)+1;
            end
        end
        if sum((arrived-delivered)==0)==M %Idle
            node_index=M+1;                                    
        elseif set_Policy==1 %The LINEAR Max Weight Policy of AoI with buffers is used
                [tilde1,node_index]=max((A(1:M)./optMU_buff).*(h-z).*((arrived-delivered)>=1));
                %[tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*(h-z).*(1-vector_successful_tx));
        elseif set_Policy==0 %The Throughput MaxWeight is used
            %if sum(max(throughput_debt(1:M),0))>0
            [tilde1,node_index]=max(A(1:M).*(arrived-delivered));
            %else
            %    node_index=mod(node_index,M)+1; %Do something similar to Round Robin
            %end        
        else
            disp('There is a problem in set_Policy. set_Whittle=1 for Max Weight, set_Whittle=0 for Greedy on P*Debt.')
            break;
        end       
    
        AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h(:);       
        h(:)=h(:)+1;                
        
        if (node_index~=(M+1)) && (rand<p(node_index)) && ((arrived(node_index)-delivered(node_index))>=1)
            delivered(node_index)=delivered(node_index)+1;
            %vector_successful_tx(node_index)=1;                
            h(node_index)=slot_index-delay_buffer(node_index,delivered(node_index))+1;
        end  
    end
    real_q(:,iteration)=delivered(1:M)/K;
end
AoI_average_GRE=mean(AoI,2); %Averaging on the ITERATIONS
optimal_cost_clients=T*AoI_average_GRE(:)/K;
optimal_cost=sum(optimal_cost_clients);
actual_throughput=mean(real_q,2);
%optimal_cost_GRE=sum(q_vector_GRE(:,1).*alpha(1:M));

end