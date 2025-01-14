function [Q,optimal_cost]=MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy,Achievable_AoI)
% function [optimal_cost_clients,optimal_cost]=MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy,Achievable_AoI)

%[actual_throughput,optimal_cost]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,work_conserving_policy)
%set_Policy=1 for LINEAR Max_Weight, set_Policy=0 for Greedy on h, set_Policy=2 for QUADRATIC Max_Weight and set_Policy=3 for LINEAR Max_Weight with NO BUFFER

%Author: Igor Kadota
%Advisor: Eytan Modiano
%Period: Fall 2018

global probability cost alpha arrival optMU optMU_NO_buff L

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
Q1=0;
V=0.1;
%% Simulation

real_q=zeros(M,num_iterations);
AoI=zeros(M,num_iterations);
node_index=0;
% Achievable_AoI=(sum(sqrt(A./p.*(3*L.^2-L)./(2.*L))))^2/M+sum(A.*3/2)/M;
for iteration=1:1:num_iterations

    %transmissions=zeros(M+1,1);
    delivered=zeros(M+1,1);    
    throughput_debt=zeros(M,1);
    h=h1*ones(M,1);
    z=z1*ones(M,1);
    Q=Q1*ones(M,1);
    wait = zeros(M,1);
    vector_successful_tx=zeros(M,1);
    QL=zeros(M,max(L));
    for slot_index=1:1:K       
                
        for node_index=1:1:M
            if vector_successful_tx(node_index)==L(node_index)
                wait(node_index)=wait(node_index)+1;
            end
            if rand<Arr(node_index) && vector_successful_tx(node_index)==0 
                % z(node_index)=0;
                vector_successful_tx(node_index)=L(node_index);% new arrival
                wait(node_index)=0;
            % else 
                % z(node_index)=z(node_index)+1;
            end

        end
        
        if sum(vector_successful_tx)==0 %Idle
            node_index=M+1;                                    
        elseif set_Policy==1 %The LINEAR Max Weight Policy is used                
                [tilde1,node_index]=max(rand(M,1).*optMU);
                % [tilde1,node_index]=max((A(1:M).*(1./optMU+p(1:M)./Arr(1:M)-p(1:M))).*(h-z).*double(vector_successful_tx>0));
                %[tilde1,node_index]=max((A(1:M)./optMU).*(h-z).*(1-vector_successful_tx));
                %[tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*(h-z).*(1-vector_successful_tx));
        elseif set_Policy==0 %The Greedy on the h is used
            %if sum(max(throughput_debt(1:M),0))>0 
            % [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*(Q+z+vector_successful_tx+1).*V).^2.*double(vector_successful_tx>0));
            % [tilde1,node_index]=max(((p(1:M).*(h-z-vector_successful_tx-1).*max(Q+z+vector_successful_tx+1-Achievable_AoI,0))).*double(vector_successful_tx>0));
            weight = zeros(M,1);
            for i=1:1:M
                weight(i)= max(QL(i,vector_successful_tx(i))+h(i)+1-Achievable_AoI(i),0).^2-max(QL(i,vector_successful_tx(i))+z(i)+vector_successful_tx(i)-Achievable_AoI(i),0).^2;
            end
            [tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*weight);
            % [tilde1,node_index]=max(p(1:M).*(h-z-vector_successful_tx-1).*double(vector_successful_tx>0));
            
            % [tilde1,node_index]=max(((V.*(2.*(z+vector_successful_tx)+1).*double(vector_successful_tx>1+((z+vector_successful_tx+1).^2-L.^2).*double(vector_successful_tx==1)))).*double(vector_successful_tx>0));
            % % [tilde1,node_index]=max(Q.*double(vector_successful_tx>0));
            %else
            %    node_index=mod(node_index,M)+1; %Do something similar to Round Robin
            %end
        elseif set_Policy==2 %The QUADRATIC Max Weight Policy is used                
            % [tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*(h+V*Q).*double(vector_successful_tx>0));  
            % [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*(h)+sqrt(A(1:M).*p(1:M)).*V.*Q).*double(vector_successful_tx>0)); 
            % [tilde1,node_index]=max(((sqrt(A(1:M).*p(1:M)).*(h-z)).^2.*(double(vector_successful_tx==1))+(V.*sqrt(A(1:M).*p(1:M)).*(Q+h+z+vector_successful_tx+1))).^2.*double(vector_successful_tx>0)); 
            % [tilde1,node_index]=max(((sqrt(A(1:M).*p(1:M)).*(h-z)).^2.*(double(vector_successful_tx==1))+2*(sqrt(A(1:M).*p(1:M)).*((h-z-vector_successful_tx-1).*max(Q+z+vector_successful_tx+1-Achievable_AoI,0))+Q.^2)).*double(vector_successful_tx>0)); 
            % [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*(h-z-1).*double(vector_successful_tx==1)+V*(Achievable_AoI-z-vector_successful_tx-1)).*double(vector_successful_tx>0)); 
            % [tilde1,node_index]=max(((sqrt(A(1:M).*p(1:M)).*h) ...
            %     +(V.*(2.*(z+vector_successful_tx)+1).*double(vector_successful_tx>1)+(z+vector_successful_tx+1).^2-L.^2).*double(vector_successful_tx==1)).*double(vector_successful_tx>0)); 
            % [tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*h.*double(vector_successful_tx>0));
            % [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*(h-z).*double(vector_successful_tx==L) ...
            %     +3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*(z+vector_successful_tx./p(1:M))./p(1:M).*double(vector_successful_tx<L) ...
            %     ).*double(vector_successful_tx>0));
%% h-z + z+L
            % [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*(h-z).*double(vector_successful_tx==L) ...
            %     +sqrt(A(1:M).*p(1:M)).*((h-z).^2-(z+1).^2).*double(vector_successful_tx==1)...
            %     +3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*(z+vector_successful_tx./p(1:M))./p(1:M).*double(vector_successful_tx<L).*double(vector_successful_tx>1) ...
            %     +3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*((z+1).^2-L.^2./p(1:M)./p(1:M)).*double(vector_successful_tx==1)).*double(vector_successful_tx>0));
            % [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*(h-z).*double(vector_successful_tx==L) ...
            %     +sqrt(A(1:M).*p(1:M)).*((h-z).^2-(z+1).^2).*double(vector_successful_tx==1)...
            %     +3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*(z+vector_successful_tx).*double(vector_successful_tx<L).*double(vector_successful_tx>1) ...
            %     +3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*((z+2).^2-L.^2).*double(vector_successful_tx==1)).*double(vector_successful_tx>0));
% \delta(O(t))=E[L(t+1)-L(t)]----->L? UPDATE of NEW PACKET?
            [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*(h-z).*double(vector_successful_tx==L) ...
                +sqrt(A(1:M).*p(1:M)).*((h-z).^2-(z+1).^2).*double(vector_successful_tx==1)...
                +3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*(z+vector_successful_tx./p(1:M))./p(1:M).*double(vector_successful_tx<L) ...
                +2.*(L-1).*sqrt(A(1:M).*p(1:M)).*(h-z)./p(1:M).*double(vector_successful_tx<L) ...
                +0*3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*((z+2).^2-L.^2./p(1:M)./p(1:M)).*double(vector_successful_tx==1)).*double(vector_successful_tx>0));



%% W+S+WS
            % [tilde1,node_index]=max((sqrt(A(1:M).*p(1:M)).*wait.*double(vector_successful_tx==L) ...
            %     +3.*(L-1).^2.*sqrt(A(1:M).*p(1:M)).*(z+vector_successful_tx./p(1:M))./p(1:M).*double(vector_successful_tx<L) ...
            %     +2.*(L-1).*sqrt(A(1:M).*p(1:M)).*wait.*double(vector_successful_tx<L)).*double(vector_successful_tx>0));



            % % [tilde1,node_index]=max(((sqrt(A(1:M).*p(1:M)).*h)+(2*(sqrt(A(1:M).*p(1:M)).*((h-z-vector_successful_tx-1).*max(Q+z+vector_successful_tx+1-Achievable_AoI,0))+Q.^2))).*double(vector_successful_tx>0)); 
            %[tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*((h+1).^2-(z+1).^2).*(1-vector_successful_tx));
        elseif set_Policy==3 %The Linear Max Weight (WITH NO BUFFER) is used                
            % [tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*(h.^2+V*Q.^2).*double(vector_successful_tx>0)); 
            [tilde1,node_index]=max(h.*double(vector_successful_tx>0));         
            %[tilde1,node_index]=max(sqrt(A(1:M).*p(1:M)).*((h+1).^2-(z+1).^2).*(1-vector_successful_tx));
        else
            disp('There is a problem in set_Policy. set_Whittle=1 for Max Weight, set_Whittle=0 for Greedy on P*Debt.')
            break;
        end       
    
        AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h(:);       
        h(:)=h(:)+1;                
        for node_index_Q = 1:1:M
            if node_index_Q ~= node_index
                Q(node_index_Q) = max(Q(node_index_Q)+h(node_index_Q)-Achievable_AoI(node_index_Q),0);
                QL(node_index_Q,vector_successful_tx(node_index_Q))= max(QL(node_index_Q,vector_successful_tx(node_index_Q))+(z(node_index_Q)+vector_successful_tx(node_index_Q)+1)-Achievable_AoI(node_index_Q),0);
            end
            % Q(node_index_Q) = Q(node_index_Q)+L(node_index);
        end
        if (node_index~=(M+1)) && (rand<p(node_index)) %&& (vector_successful_tx(node_index)>=1)
            vector_successful_tx(node_index)=vector_successful_tx(node_index)-1;
            Q(node_index) = max(Q(node_index)+(h(node_index)+1)-Achievable_AoI(node_index),0);
            if vector_successful_tx(node_index)>0
                QL(node_index,vector_successful_tx(node_index))= max(QL(node_index,vector_successful_tx(node_index))+(z(node_index)+vector_successful_tx(node_index)+1)-Achievable_AoI(node_index),0);            
            else
                h(node_index)=z(node_index)+1;
                z(node_index)=0;
                % Q(node_index)=0;
                delivered(node_index)=delivered(node_index)+1;
            end
        else
            Q(node_index) = max(Q(node_index)+h(node_index)-Achievable_AoI(node_index),0);
            QL(node_index,vector_successful_tx(node_index))= max(QL(node_index,vector_successful_tx(node_index))+(z(node_index)+vector_successful_tx(node_index)+1)-Achievable_AoI(node_index),0);            
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
actual_throughput=mean(real_q,2);
%optimal_cost_GRE=sum(q_vector_GRE(:,1).*alpha(1:M));

end