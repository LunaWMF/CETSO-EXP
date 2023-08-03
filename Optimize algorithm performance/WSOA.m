% War Strategy Optimization Algorithm
function [King_fit,King,Convergence_curve]=WSOA(Soldiers_no,Max_iter,lb,ub,dim,fobj)

King=zeros(1,dim);
King_fit=inf;
Positions=initialization(Soldiers_no,dim,ub,lb);
pop_size=size(Positions,1);Convergence_curve=zeros(1,Max_iter);

Positions_new=zeros(size(Positions));
fitness_old=inf*ones(1,pop_size);
fitness_new=inf*ones(1,pop_size);
l=1;
W1=2*ones(1,pop_size);
Wg=zeros(1,pop_size);
Trajectories=zeros(Soldiers_no,Max_iter);
R=0.8;
for  j=1:size(Positions,1)
         fitness=fobj(Positions(j,:));
         fitness_old(j)=fitness;

        if fitness<King_fit
            King_fit=fitness;
            King=Positions(j,:);
        end
        
end
 
while l<Max_iter
[~,tindex]=sort(fitness_old);
Co=Positions(tindex(2),:);
  
  iter =l;
   com=randperm(pop_size);  
    for i=1:pop_size
        RR=rand;
        if RR<R

             D_V(i,:)=2*RR*(King-Positions(com(i),:))+1*W1(i)*rand*(Co-Positions(i,:)); 
        else
            D_V(i,:)=2*RR*(Co-King)+1*rand*(W1(i)*King-Positions(i,:)); 
        end
            Positions_new(i,:)=Positions(i,:)+D_V(i,:);
            Flag4ub=Positions_new(i,:)>ub;
            Flag4lb=Positions_new(i,:)<lb;
            Positions_new(i,:)=(Positions_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
            fitness=fobj(Positions_new(i,:));
            fitness_new(i) = fitness;
        if fitness<King_fit
            King_fit=fitness;
            King=Positions_new(i,:);
        end
        if fitness<fitness_old(i)
            Positions(i,:)=Positions_new(i,:);
            fitness_old(i)=fitness;
            Wg(i)=Wg(i)+1;
            W1(i)=1*W1(i)*(1-Wg(i)/Max_iter)^2;
        end
    end
    if l<1000
        [~,tindex1]=max(fitness_old);
        Positions(tindex1,:)=lb+rand*(ub-lb);
    end
    l=l+1;   
    Convergence_curve(l)=King_fit;
end



