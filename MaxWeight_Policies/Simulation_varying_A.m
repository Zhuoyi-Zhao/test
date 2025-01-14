function [K,vector_A,num_iterations,set_simplified,optimal_cost_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Whittle,optimal_cost_Rand,optimal_cost_Rand_FIFO,optimal_cost_Whittle_NO_Buffer,optimal_cost_Max_Weight_FIFO,optimal_cost_Rand_NO_Buffer,optimal_cost_Throughput_Max_Weight_FIFO,optimal_cost_Max_Weight_NO_buff]=Simulation_varying_A
%[K,vector_A,num_iterations,set_simplified,optimal_cost_Linear_Max_Weight,optimal_cost_Buffered_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Throughput_Max_Weight,optimal_cost_Rand]=Simulation_varying_A
%[K,vector_A,num_iterations,set_simplified,optimal_cost_Linear_Max_Weight,optimal_cost_Quadratic_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Whittle_Single_Buffer,optimal_cost_Rand,optimal_cost_Whittle_NO_Buffer,throughput_Linear_Max_Weight,throughput_Quadratic_Max_Weight,throughput_Whittle_Single_Buffer,throughput_Rand,throughput_Whittle_NO_Buffer]=Simulation_varying_A
%[K,vector_M,num_iterations,set_simplified,optimal_cost_Linear_Max_Weight,optimal_cost_Quadratic_Max_Weight,optimal_cost_LowerBound_stoch,optimal_cost_Whittle_Single_Buffer_Throughput_Constraint,average_throughput_Linear_Max_Weight,average_throughput_Quadratic_Max_Weight,average_throughput_Whittle_Single_Buffer]=Simulation_varying_A

clear
clc

global probability cost alpha arrival optQ optMU optMU_buff optMU_NO_buff L

tic;

%% Model Setup
M=4;
vector_A=0.01:0.01:0.4;
%epsilon=0.9;
%vector_M=[2,4,6];
L = [5;5;5;5];
num_setups=length(vector_A);
%num_iterations=10;
num_iterations=10;
T=1;
set_simplified=1;
h1=1; %initial age
z1=0; %HAS TO BE LOWER THAN H1

%% Parameters Setup
optimal_cost_Rand=zeros(1,num_setups);
% optimal_cost_Rand_FIFO=zeros(1,num_setups);
% optimal_cost_Rand_NO_Buffer=zeros(1,num_setups);
% optimal_cost_Whittle_NO_Buffer=zeros(1,num_setups);
% optimal_cost_Whittle=zeros(1,num_setups);
% optimal_cost_Max_Weight_FIFO=zeros(1,num_setups);
% optimal_cost_Max_Weight_NO_buff=zeros(1,num_setups);
optimal_cost_Max_Weight=zeros(1,num_setups);
% optimal_cost_Throughput_Max_Weight_FIFO=zeros(1,num_setups);
optimal_cost_LowerBound_stoch=zeros(1,num_setups);

%throughput_Linear_Max_Weight=zeros(max(vector_A),num_setups);
%throughput_Quadratic_Max_Weight=zeros(max(vector_A),num_setups);
%throughput_Whittle_Single_Buffer=zeros(max(vector_A),num_setups);
%throughput_Rand=zeros(max(vector_A),num_setups);
%throughput_Whittle_NO_Buffer=zeros(max(vector_M),num_setups);
%throughput_Buffered_Max_Weight=zeros(max(vector_M),num_setups);
%throughput_Throughput_Max_Weight=zeros(max(vector_M),num_setups);


%% Simulations
for count=1:num_setups
    count    
    %K=1000000*M;
    K=500*M; 
    %K=500*M;
    cost=zeros(M,1);
    alpha=[4;4;1;1];
    arrival=vector_A(count)*[1;3/4;1/2;1/4];
    probability=[1/4;1/2;3/4;1];
    
    p=zeros(M,1);
    p=p+probability; %Making sure that "probability" has the right size

    A=zeros(M,1);
    A =A+alpha; %Making sure that "alpha" has the right size
    
    Arr=zeros(M,1);
    Arr=Arr+arrival; %Making sure that "constraint" has the right size
    
    Lower_Bound_find_q(M); %findind the optQ            
    optimal_cost_LowerBound_stoch(count)=sum(alpha.*(1+1./optQ))/(2*M);
    disp('Lower Bound for Stochastic Arrivals')         
    
    optMU=sqrt(alpha./probability);
    sum_optMU=sum(optMU);
    optMU=optMU/sum_optMU;
        
    optMU_NO_buff=sqrt(alpha./(arrival.*probability));
    sum_optMU_NO_buff=sum(optMU_NO_buff);
    optMU_NO_buff=optMU_NO_buff/sum_optMU_NO_buff;
    
    Random_find_mu_buffered_BerBer1(M) %findind the optMU_buff   
    optimal_cost_Rand(count)=(sum(sqrt(A./p.*(2*L.^2+L-1)./(2.*L))))^2/M+(sum(A.*(1./Arr-1)))/M;
    
    set_Policy=1; %LINEAR MaxWeight    
    [aux,optimal_cost_Max_Weight(count)]=MaxWeight_and_Greedy_simulation_stochastic(K,T,M,h1,z1,num_iterations,set_Policy);    
    optimal_cost_Max_Weight(count)=optimal_cost_Max_Weight(count)/(M*T);
    %average_throughput_Linear_Max_Weight(count)=mean(throughput_Linear_Max_Weight);    
    disp('Linear Max-Weight policy with Optimal Buffer')
    
    toc
end

figure(3)
hold on
[Rand_plot]=plot(vector_A,optimal_cost_Rand,'bo--','LineWidth',3,'MarkerSize',20);
[Linear_Max_Weight_plot]=plot(vector_A,optimal_cost_Max_Weight,'rx--','LineWidth',3,'MarkerSize',25);
[Lower_Bound]=plot(vector_A,optimal_cost_LowerBound_stoch,'k-','LineWidth',3,'MarkerSize',20);
legend([Rand_plot,Linear_Max_Weight_plot,Lower_Bound],'Random. with Single pkt queues','Max-Weight with Single pkt queues''Lower Bound','Location','NorthWest');
ylabel('Expected Weighted Sum AoI')
xlabel('Arrival Rate, \lambda')
%title('Comparing the performance of different policies')
hold off

toc;
end