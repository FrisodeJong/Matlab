%randomize lambda which stays fixed after that
rho=rand(20,1);
Lambda=zeros(20,20);
for i=1:20
    for j=1:20
        Lambda(i,j)=2*rho(i)*(1-rho(j));
    end
end
%randomize D demand which changes each iteration
D=zeros(20,20);
for i=1:20
    for j=1:20
        D(i,j)=poissrnd(Lambda(i,j));
    end
end
%randomize distances
a = 100;
b = 1500;
M = (b-a).*rand(20,20) + a;
for i=1:20
    M(i,i)=0;
end
%rewards
R=zeros(20,20);
for i=1:20
    for j=1:20
        R(i,j)=rho(i)*M(i,j);
    end
end
%costs
C=zeros(20,20);
for i=1:20
    for j=1:20
        C(i,j)=1.2*M(i,j);
    end
end
%initialize while loop
n=1;
N=1000;
%initinalize starting value function approximation, alpha approximation
%V=1+zeros(20,1);
V=averageV;
alpha=10/(9+n);
cityhistory=zeros(22,1);
safecityhistory=zeros(22,1000);
safeV=zeros(20,1000);
while n<N
    %V=1+zeros(20,1);
    V=averageV;
    %initialize starting city
    S=[1:1:20];
    city=S(1,1);     
    %choose random sample path
    D=zeros(20,20);
    for i=1:20
        for j=1:20
            D(i,j)=poissrnd(Lambda(i,j));
        end
    end 
    u=1;
    for t=0:21
        cityhistory(u,1)=city;
        u=u+1;
        %solve for x based on random demands    
        D=zeros(20,20);
        for i=1:20
            for j=1:20
                D(i,j)=poissrnd(Lambda(i,j));
            end
        end 
        i=city;  
        f1=zeros(20,1);
        f2=zeros(20,1);
        for p=1:20
            f1(p,1)=-R(i,p)-V(p,1); 
        end
        for q=1:20
            f2(q,1)=C(i,q)-V(q,1);
        end
        f=[f1;f2];
        intcon=[1:1:40];
        A=eye(40,40);
        b1=D(i,:).';
        b2=ones(20,1);
        b=[b1;b2]; %constraint of x^L_ij being 0 whenever d_ij = 0
        Aeq=ones(1,40);
        beq=1; %making sure sum of all x's equals 1
        lb=zeros(40,1);
        ub=ones(40,1);
        x=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
        %solve for v given calculated optimal action
        f1=zeros(20,1);
        f2=zeros(20,1);
        for p=1:20
            f1(p,1)=R(i,p)+V(p,1); 
        end
        for q=1:20
            f2(q,1)=-C(i,q)+V(q,1);
        end
        f=[f1;f2];
        v=f.'*x;
        %update V in current city      
        V(i,1) = (1-alpha)*V(i,1) + alpha*v; 
        
        %compute where trucker will be in the next time epoch
        for i=1:40  
            if x(i,1)~=0
                if i<=20
                    city=i;
                else
                    city=i-20;
                end
            end
        end
        
    end
    safecityhistory(:,n)=cityhistory;
    safeV(:,n)=V;
    n=n+1;
    alpha=10/(9+n);  
end
%average Value after 1000 iterations
for i=1:20
    averageV(i,1)=mean(safeV(i,:));
end
%city counter
for i=1:20
    Citycounter(i,1)=sum(safecityhistory(:)==i);
end


