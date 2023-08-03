%  Sine Cosine Algorithm (SCA)2016

function [Destination_fitness,Destination_position,Convergence_curve]=SCA(N,Max_iteration,lb,ub,dim,fobj)
X=initialization(N,dim,ub,lb);
Destination_position=zeros(1,dim);
Destination_fitness=inf;
Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));

for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end
t=2; 
while t<=Max_iteration
    a = 2;
    Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); 
    for i=1:size(X,1) 
        for j=1:size(X,2) 
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            if r4<0.5
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
        end
    end
    
    for i=1:size(X,1)
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        Objective_values(1,i)=fobj(X(i,:));
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    Convergence_curve(t)=Destination_fitness;
    t=t+1;
end
end
%%
function X=initialization(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); 
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end