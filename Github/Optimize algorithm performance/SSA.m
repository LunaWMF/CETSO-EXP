% (Salp Swarm Algorithm: SSA))2017
function [FoodFitness,FoodPosition,Convergence_curve]=SSA(N,Max_iter,lb,ub,dim,fobj)
SalpPositions=initialization(N,dim,ub,lb);
if size(ub,2)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
else
    ub=ub';
    lb=lb';
end
Convergence_curve = zeros(1,Max_iter);
FoodPosition=zeros(1,dim);
FoodFitness=inf;
for i=1:size(SalpPositions,1)
    SalpFitness(1,i)=fobj(SalpPositions(i,:));
end
[sorted_salps_fitness,sorted_indexes]=sort(SalpFitness);
for newindex=1:N
    Sorted_salps(newindex,:)=SalpPositions(sorted_indexes(newindex),:);
end

FoodPosition=Sorted_salps(1,:);
FoodFitness=sorted_salps_fitness(1);
l=2; 
while l<Max_iter+1
    c1 = 2*exp(-(4*l/Max_iter)^2);
    for i=1:size(SalpPositions,1)
        SalpPositions= SalpPositions';
        if i<=N/2
            for j=1:1:dim
                c2=rand();
                c3=rand();
                if c3<0.5 
                    SalpPositions(j,i)=FoodPosition(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                else
                    SalpPositions(j,i)=FoodPosition(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                end
            end
        elseif i>N/2 && i<N+1
            point1=SalpPositions(:,i-1);
            point2=SalpPositions(:,i);
            SalpPositions(:,i)=(point2+point1)/2;
        end
        SalpPositions= SalpPositions';
    end
    for i=1:size(SalpPositions,1)
        Tp=SalpPositions(i,:)>ub';Tm=SalpPositions(i,:)<lb';
        SalpPositions(i,:)=(SalpPositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
        SalpFitness(1,i)=fobj(SalpPositions(i,:));
        if SalpFitness(1,i)<FoodFitness
            FoodPosition=SalpPositions(i,:);
            FoodFitness=SalpFitness(1,i);
        end
    end
    Convergence_curve(l)=FoodFitness;
    l = l + 1;
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
