function [actual_throughput,optimal_cost]=Whittle_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy)
%[actual_throughput,optimal_cost]=Whittle_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy)
%set_Policy=1 for NO BUFFER (Yu-Pin) case. set_Policy=0 for "KEEP THE LAST PACKET" BUFFER case.

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
    %throughput_debt=zeros(M,1);
    h=h1*ones(M,1);
    z=z1*ones(M,1);
    vector_successful_tx=zeros(M,1);

    for slot_index=1:1:K       
                
        for node_index=1:1:M
            if rand<Arr(node_index)
                z(node_index)=0;
                vector_successful_tx(node_index)=0;% new arrival
            else
                z(node_index)=z(node_index)+1;
            end
        end
        if sum(vector_successful_tx)==M %Idle
            node_index=M+1;
        elseif set_Policy==1 %NO BUFFER (Yu-Pin) case.            
            if sum(z==0)==0
                node_index=M+1;
            else
                [tilde1,node_index]=max((h.^2./2+(1./Arr(1:M)-1/2).*h).*(z==0));                                         
            end
        elseif set_Policy==0 %"KEEP THE LAST PACKET" BUFFER case.            
            %x=zeros(M,1);
            x=((h-z+1)+Arr(1:M).*(z.*(z+1))./2)./(1-Arr(1:M)+(z+1).*Arr(1:M));
            wittle_index=zeros(M,1);
            for i=1:M
                if h(i)<=(Arr(i)*(z(i)+1)^2/2+(2-Arr(i)/2)*(z(i)+1))
                    wittle_index(i)=(h(i)-z(i)-1)/Arr(i);
                else
                    wittle_index(i)=x(i)^2/2+(1/Arr(i)-1/2)*x(i);
                end
            end
            [tilde1,node_index]=max(wittle_index.*(1-vector_successful_tx));
        else
            disp('There is a problem in set_Policy. set_Whittle=1 for Max Weight, set_Whittle=0 for Greedy on P*Debt.')
            break;
        end       
    
        AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h(:);       
        h(:)=h(:)+1;                
        
        if (node_index~=(M+1)) && (rand<p(node_index)) && (vector_successful_tx(node_index)==0)
            delivered(node_index)=delivered(node_index)+1;
            vector_successful_tx(node_index)=1;                
            h(node_index)=z(node_index)+1;
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