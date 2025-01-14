function Random_find_mu_buffered_BerBer1(M)
%This is the algrothm that find the optimal scheduling probability for
%randomized policies.

global probability alpha arrival optMU_buff

iterations=1000;

epsilon=0.000001;

p=zeros(M,1);
p=p+probability; %Making sure that "probability" has the right size

A=zeros(M,1);
A =A+alpha; %Making sure that "alpha" has the right size

%C=zeros(M,1);
%C =[C+cost]; %Making sure that "cost" has the right size

Arr=zeros(M,1);
Arr=Arr+arrival; %Making sure that "constraint" has the right size

%GammaI=(A.*p)./(M.*Q.^2);
J_prime=zeros(M,1);
mu_initial=(Arr+epsilon)./p;
for node_index=1:M
    [J_prime(node_index)]=g(mu_initial(node_index),node_index,M);
end
GammaI=-(A./M).*J_prime;
maxGamma=max(GammaI);
minGamma=0;
%S=sum(mu_initial);
for count=1:iterations    
    Gamma_aux=(maxGamma+minGamma)/2;
    mu_aux=zeros(M,1);
    mu_inv=zeros(M,1);
    for node=1:M
        mu_inv(node)=g_inv(-(Gamma_aux*M/A(node)),node,M);
        mu_aux(node)=max([mu_initial(node),mu_inv(node)]);
    end    
    S=sum(mu_aux);
    if S<=1
        maxGamma=Gamma_aux;
    else
        minGamma=Gamma_aux;
    end
end
optMU_buff=zeros(M,1);
for node=1:M
    optMU_buff(node)=max([mu_initial(node),g_inv(-(maxGamma*M/A(node)),node,M)]);                        
end    
if sum(Arr./p)>=1
    disp('Infeasible system');     
    optMU_buff=epsilon*ones(M,1);
end

end

function [J_prime]=g(u,node_index,M)

global probability arrival

p=zeros(M,1);
p=p+probability; %Making sure that "probability" has the right size

%A=zeros(M,1);
%A =A+alpha; %Making sure that "alpha" has the right size

Arr=zeros(M,1);
Arr=Arr+arrival; %Making sure that "constraint" has the right size

J_prime=(Arr(node_index)/(p(node_index)*u^2))*(2/(p(node_index)*u)-1)-(p(node_index)*(1-Arr(node_index)))/(p(node_index)*u-Arr(node_index))^2;

end 

function [mu]=g_inv(J_prime,node_index,M)

global probability alpha arrival

p=zeros(M,1);
p=p+probability; %Making sure that "probability" has the right size

A=zeros(M,1);
A =A+alpha; %Making sure that "alpha" has the right size

Arr=zeros(M,1);
Arr=Arr+arrival; %Making sure that "constraint" has the right size

iterations=1000;

mu_min=Arr(node_index)/p(node_index);
mu_max=1;

for count=1:iterations
    mu_aux=(mu_max+mu_min)/2;    
    J_prime2=(Arr(node_index)/(p(node_index)*mu_aux^2))*(2/(p(node_index)*mu_aux)-1)-(p(node_index)*(1-Arr(node_index)))/(p(node_index)*mu_aux-Arr(node_index))^2;
    result=J_prime-J_prime2;    
    if result<0
        mu_max=mu_aux;
    else
        mu_min=mu_aux;
    end
end
mu=mu_min;
end