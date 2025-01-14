% function [K,vector_p,num_iterations,set_simplified,optimal_cost_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Whittle,optimal_cost_Rand,optimal_cost_Rand_FIFO,optimal_cost_Whittle_NO_Buffer,optimal_cost_Max_Weight_FIFO,optimal_cost_Rand_NO_Buffer,optimal_cost_Throughput_Max_Weight_FIFO,optimal_cost_Max_Weight_NO_buff]=Simulation_varying_p
%[K,vector_A,num_iterations,set_simplified,optimal_cost_Linear_Max_Weight,optimal_cost_Buffered_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Throughput_Max_Weight,optimal_cost_Rand]=Simulation_varying_A
%[K,vector_A,num_iterations,set_simplified,optimal_cost_Linear_Max_Weight,optimal_cost_Quadratic_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Whittle_Single_Buffer,optimal_cost_Rand,optimal_cost_Whittle_NO_Buffer,throughput_Linear_Max_Weight,throughput_Quadratic_Max_Weight,throughput_Whittle_Single_Buffer,throughput_Rand,throughput_Whittle_NO_Buffer]=Simulation_varying_A
%[K,vector_M,num_iterations,set_simplified,optimal_cost_Linear_Max_Weight,optimal_cost_Quadratic_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Whittle_Single_Buffer_Throughput_Constraint,average_throughput_Linear_Max_Weight,average_throughput_Quadratic_Max_Weight,average_throughput_Whittle_Single_Buffer]=Simulation_varying_A

clear
clc

global probability cost alpha arrival optQ optMU optMU_buff optMU_NO_buff L

tic;

%% Model Setup
M=2;
vector_p=0.2:0.05:0.5;
%epsilon=0.9;
%vector_M=[2,4,6];
L = [3;30];
num_setups=length(vector_p);
num_iterations=10;
%num_iterations=2;
T=1;
set_simplified=1;
h1=1; %initial age
z1=0; %HAS TO BE LOWER THAN H1

%% Parameters Setup
optimal_cost_Rand=zeros(1,num_setups);
optimal_cost_Rand_NoSwitch = zeros(1,num_setups);
optimal_cost_Rand_NoSwitch_Simulated = zeros(1,num_setups);
optimal_cost_Max_Weight_Greedy=zeros(1,num_setups);
optimal_cost_Max_Weight_LinearDPP=zeros(1,num_setups);
optimal_cost_Max_Weight_Age_debt=zeros(1,num_setups);
optimal_cost_LowerBound_stoch=zeros(1,num_setups);



%% Simulations
for count=1:num_setups
    count    
    %K=1000000*M;
    K=50000*M; 
    %K=500*M;
    
    cost=zeros(M,1);
    alpha=[5;1];
    arrival=[1;1];
    probability=vector_p(count)*[1;1];
    % probability=[vector_p(count);0.4;0.6;0.8]; 
    p=zeros(M,1);
    p=p+probability; %Making sure that "probability" has the right size

    A=zeros(M,1);
    A =A+alpha; %Making sure that "alpha" has the right size
    
    Arr=zeros(M,1);
    Arr=Arr+arrival; %Making sure that "constraint" has the right size
    
    Lower_Bound_find_q(M); %findind the optQ            
    optimal_cost_LowerBound_stoch(count)=sum(alpha.*(1+L./optQ))/(2*M);
    disp('Lower Bound for Stochastic Arrivals')         
    
    optMU=sqrt(A./p.*(3*L.^2-L)./(2.*L));
    sum_optMU=sum(optMU);
    optMU=optMU/sum_optMU;
        
    optMU_NO_buff=sqrt(alpha./(arrival.*probability));
    sum_optMU_NO_buff=sum(optMU_NO_buff);
    optMU_NO_buff=optMU_NO_buff/sum_optMU_NO_buff;
    
    E_Y2 = zeros(M,1);
    for j=1:M
        E_Y2(j) = 1 + 2*(L(j)-1) + ((L(j)-1)*(L(j)-p(j)))/p(j);
    end

    Random_find_mu_buffered_BerBer1(M) %findind the optMU_buff   
    optimal_cost_Rand(count)=(sum(sqrt(A./p.*(3*L.^2-L)./(2.*L))))^2/M+sum(A.*3/2)/M;
    % No Switching Randomized policy
    mu_opt_no_switch = optimize_mu_for_random_policy(M, p, L, alpha, E_Y2);
    optimal_cost_Rand_NoSwitch(count) = weighted_AoI(mu_opt_no_switch, p, L, E_Y2, alpha)-1.5
    [~,optimal_cost_Rand_NoSwitch_Simulated(count)] = No_Switch_Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,mu_opt_no_switch)

    % set_Policy=3; %1: Randomized 0: Age-Debt,Q 2: DPP,h+Q  3:QuaDPP, h.^2+Q.^2  
    % % [aux,optimal_cost_Max_Weight(count)]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy);    
    % Achievable_AoI = (1./p./optMU.*(3*L.^2-L)./(2.*L))./1+1.*3/2;
    % [optcostgreedy,optimal_cost_Max_Weight_Greedy(count)]=MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy,Achievable_AoI);    
    % optimal_cost_Max_Weight_Greedy(count)=optimal_cost_Max_Weight_Greedy(count)/(M*T)
    % %average_throughput_Linear_Max_Weight(count)=mean(throughput_Linear_Max_Weight);    
    % disp('Linear Max-Weight policy with Optimal Buffer')
    % toc 
    % 
    % % set_Policy=0; %1: Randomized 0: Age-Debt,Q 2: DPP,h+Q  3:QuaDPP, h.^2+Q.^2  
    % % % [aux,optimal_cost_Max_Weight(count)]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy);    
    % % Achievable_AoI = (1./p./optMU.*(3*L.^2-L)./(2.*L))./2.2+1.*3/2;
    % % % Achievable_AoI = (1 ./ p ./ optMU .* (3 * L.^2 - L) ./ (2 .* L))/2 + 1 * 3 / 2;
    % % [aux,optimal_cost_Max_Weight_Age_debt(count)]=MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy,Achievable_AoI);    
    % % optimal_cost_Max_Weight_Age_debt(count)=optimal_cost_Max_Weight_Age_debt(count)/(M*T)
    % % %average_throughput_Linear_Max_Weight(count)=mean(throughput_Linear_Max_Weight);    
    % % disp('Linear Max-Weight policy with Optimal Buffer')
    % 
    % 
    % set_Policy=2; %1: Randomized 0: Age-Debt,Q 2: DPP,h+Q  3:QuaDPP, h.^2+Q.^2  
    % % [aux,optimal_cost_Max_Weight(count)]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy);    
    % Achievable_AoI = (1./p./optMU.*(3*L.^2-L)./(2.*L))./1+1.*3/2;
    % [optcostgreedy,optimal_cost_Max_Weight_QuaDPP(count)]=MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy,Achievable_AoI);    
    % optimal_cost_Max_Weight_QuaDPP(count)=optimal_cost_Max_Weight_QuaDPP(count)/(M*T)
    % %average_throughput_Linear_Max_Weight(count)=mean(throughput_Linear_Max_Weight);    
    % disp('Linear Max-Weight policy with Optimal Buffer')
    % toc  
    % 
    % set_Policy=0; %1: Randomized 0: Age-Debt,Q 2: DPP,h+Q  3:QuaDPP, h.^2+Q.^2  
    % % [aux,optimal_cost_Max_Weight(count)]=Random_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy);    
    % Achievable_AoI = (1./p./optMU.*(3*L.^2-L)./(2.*L))./3+1.*3/2;
    % [aux,optimal_cost_Max_Weight_Age_debt(count)]=LQ_MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy,Achievable_AoI);    
    % optimal_cost_Max_Weight_Age_debt(count)=optimal_cost_Max_Weight_Age_debt(count)/(M*T)
    % %average_throughput_Linear_Max_Weight(count)=mean(throughput_Linear_Max_Weight);    
    % disp('Linear Max-Weight policy with Optimal Buffer')
    % toc
end

figure(4)
hold on
[Rand_NoSwitch_plot] = plot(vector_p,optimal_cost_Rand_NoSwitch,'cx-','LineWidth',3,'MarkerSize',15);
[Rand_NoSwitch_Simulated_plot] = plot(vector_p,optimal_cost_Rand_NoSwitch_Simulated,'bx','LineWidth',3,'MarkerSize',15);
[Rand_plot]=plot(vector_p,optimal_cost_Rand,'ro-','LineWidth',3,'MarkerSize',20);
% [Linear_Max_Weight_plot_Age_debt]=plot(vector_p,optimal_cost_Max_Weight_Age_debt,'mx--','LineWidth',3,'MarkerSize',25);
% [Linear_Max_Weight_plot_Greedy]=plot(vector_p,optimal_cost_Max_Weight_Greedy,'g>--','LineWidth',3,'MarkerSize',25);
% [Linear_Max_Weight_plot_QuaDPP]=plot(vector_p,optimal_cost_Max_Weight_QuaDPP,'r^--','LineWidth',3,'MarkerSize',25);
[Lower_Bound]=plot(vector_p,optimal_cost_LowerBound_stoch,'k-','LineWidth',3,'MarkerSize',20);
% legend([Rand_plot, Rand_NoSwitch_Simulated_plot,Linear_Max_Weight_plot_Age_debt, Linear_Max_Weight_plot_Greedy, Linear_Max_Weight_plot_QuaDPP, Lower_Bound, Rand_NoSwitch_plot],...
%     'Random (Switch)','Random (No-Switch)simulated','Age Debt','Greedy','QuaDPP','Lower Bound','Random (No-Switch)','Location','NorthWest');ylabel('Expected Weighted Sum AoI')
legend([Rand_plot, Rand_NoSwitch_Simulated_plot,Rand_NoSwitch_plot,Lower_Bound],...
    'Switching Randomized Policy','Simulated No-Switching Randomized Policy','No-Switching Randomized Policy','Lower Bound','Location','NorthWest');ylabel('Expected Weighted Sum AoI')
xlabel('Channel Reliabilities')
%title('Comparing the performance of different policies')
hold off

toc;
