function Lower_Bound_find_q(M)

global probability alpha arrival optQ

iterations=200;

p=zeros(M,1);
p=p+probability; %Making sure that "probability" has the right size

A=zeros(M,1);
A =A+alpha; %Making sure that "alpha" has the right size

%C=zeros(M,1);
%C =[C+cost]; %Making sure that "cost" has the right size

Arr=zeros(M,1);
Arr =Arr+arrival; %Making sure that "arrival" has the right size

optQ=zeros(M,1);

if sum(Arr./p)<=1
    for node=1:M
        optQ(node)=Arr(node);
    end
else    
    Gamma_star=(sum(sqrt(A./p)))^2/M;
    QI=sqrt(A.*p./(M*Gamma_star));
    if sum((QI-Arr)>0)==0
        for node=1:M
            optQ(node)=QI(node);
        end
    else
        GammaI=A.*p./(M*(Arr.^2));
        maxGamma=max(GammaI);
        minGamma=0;
        %Q_initial=sqrt(A./(M*maxGamma));
        %S=sum(mu_initial);
        for count=1:iterations
            Gamma_aux=(maxGamma+minGamma)/2;
            Q_aux=zeros(M,1);
            for node=1:M
                Q_aux(node)=min([Arr(node),sqrt(A(node).*p(node)./(M*Gamma_aux))]);
            end    
            S=sum(Q_aux);
            if S<=1
                maxGamma=Gamma_aux;
            else
                minGamma=Gamma_aux;
            end
        end
        optQ=zeros(M,1);
        for node=1:M
            optQ(node)=min([Arr(node),sqrt(A(node).*p(node)./(M*maxGamma))]);
        end    
    end    
end
end