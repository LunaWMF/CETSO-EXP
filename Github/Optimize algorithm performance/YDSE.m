
%Young's Double-Slit Experiment (YDSE) Optimizer
function [Best_Fringe_fitness,BestFringe,cgcurve]=YDSE(NP,Max_iter,LB,UB,Dim,F_obj)

L=1; 
d=5*10^-3; 
I=0.01;
Lambda=5*10^-6; 
Delta=0.38;
S=initialization(NP,Dim,UB,LB);FS=zeros(NP,Dim);
SS=zeros(NP,Dim);
S_Mean=mean(S);
for i=1:NP
    FS(i,:)= S(i,:)+I.*(2*rand-1).*(S_Mean-S(i,:)) ; 
    SS(i,:)= S(i,:)-I.*(2*rand-1).*(S_Mean-S(i,:)) ; 
end
BestFringe=zeros(1,Dim);   
Best_Fringe_fitness = inf;  
Fringes_fitness=zeros(1,NP);
X=zeros(NP,Dim);
for i=1:NP
    if i==1 
        Delta_L=0;
        for j=1:Dim
            X(i,j)=(FS(i,j)+SS(i,j))/2+ Delta_L;
        end
    elseif mod(i,2)==0
        m=i-1;
        Delta_L=(2*m+1)*Lambda/2;
        for j=1:Dim
            X(i,j)= (FS(i,j)+SS(i,j))/2+Delta_L;
        end
    else          
        m=i-1;
        Delta_L=(m)*Lambda;
        for j=1:Dim
            X(i,j)= (FS(i,j)+SS(i,j))/2+Delta_L;
        end
    end
end
for i=1:NP
    F_UB=X(i,:)>UB;
    F_LB=X(i,:)<LB;
    X(i,:)=(X(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
end
for i=1:NP
    Fringes_fitness(1,i) = F_obj(X(i,:));
end
[Fringes_fitness,sorted_indices]=sort(Fringes_fitness);
temp_X=X;
X=temp_X(sorted_indices,:);
Best_Fringe_fitness=Fringes_fitness(1,1);
BestFringe = X(1,:);
Int_max0=10^-20;
iter=1;
X_New=zeros(NP,Dim);
while iter<Max_iter+1
    cgcurve(iter)=Best_Fringe_fitness;
    [Fringes_fitness,sorted_indices]=sort(Fringes_fitness);
    temp_X=X;
    X=temp_X(sorted_indices,:);
    q=iter/Max_iter;  
    Int_max=Int_max0.*q;  
    for  i=1:NP
        a=iter.^(2*rand()-1); 
        H=2.*rand(1,Dim)-1;
        z=a./H;  
        r1=rand;
        
        if i==1 
            beta=q.*cosh(pi./iter); 
            A_bright=2./(1+sqrt(abs(1-beta.^2))); 
            A=1:2:NP;
            randomIndex = randi(length(A));
            X_New(i,:)=BestFringe+(Int_max.*A_bright.*X(i,:)-r1.*z.*X(A(randomIndex),:)); 
        elseif  mod(i,2)==1 
            beta=q.*cosh(pi./iter);
            A_bright=2./(1+sqrt(abs(1-beta.^2))); 
            m=i-1;
            y_bright=Lambda*L*m/d; 
            Int_bright=Int_max.*cos((pi*d)/(Lambda*L)*y_bright).^2;   
            s=randperm(NP,2);
            g=2*rand-1;
            Y=(X(s(2),:)-X(s(1),:)); 
            X_New(i,:)=X(i,:)-((1-g).*A_bright.*Int_bright.*X(i,:)+g.*Y); 
        else 
            A_dark=Delta*atanh(-(q)+1); 
            m=i-1;
            y_dark=Lambda*L*(m+0.5)/d; 
            Int_dark=Int_max.*cos((pi*d)/(Lambda*L)*y_dark).^2;
            X_New(i,:)=X(i,:)-(r1.*A_dark.*Int_dark.*X(i,:)-z.*BestFringe); 
        end
        F_UB=X_New(i,:)>UB;
        F_LB=X_New(i,:)<LB;
        X_New(i,:)=(X_New(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
    end
    X_New_Fitness=zeros(1,size(X,1));
    for i=1:NP
        X_New_Fitness(1,i) = F_obj(X_New(i,:));
        if X_New_Fitness(1,i)<Fringes_fitness(1,i)
            X(i,:)=X_New(i,:);
            Fringes_fitness(1,i)=X_New_Fitness(1,i);
            if Fringes_fitness(1,i) < Best_Fringe_fitness
                Best_Fringe_fitness=Fringes_fitness(1,i);
                BestFringe = X(i,:);
            end
        end
    end
    iter=iter+1;
end
end