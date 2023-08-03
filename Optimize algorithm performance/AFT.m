%% Ali baba and the Forty Thieves (AFT) 

function [fitness,gbest,ccurve]=AFT(noThieves,itemax,lb,ub,dim,fobj)

ccurve=zeros(1,itemax);

if size(ub,1)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
end

xth=zeros(noThieves, dim);
for i=1:noThieves 
    for j=1:dim
        xth(i,j)=  lb(j)-rand()*(lb(j)-ub(j)) ;
    end
end

fit=zeros(noThieves, 1);

for i=1:noThieves
     fit(i,1)=fobj(xth(i,:));
end

fitness=fit; 
[sorted_thieves_fitness,sorted_indexes]=sort(fit);

for index=1:noThieves
    Sorted_thieves(index,:)=xth(sorted_indexes(index),:);
end

gbest=Sorted_thieves(1,:); 
fit0=sorted_thieves_fitness(1);

best = xth;        
xab=xth;        

for ite=1:itemax

    Pp=0.1*log(2.75*(ite/itemax)^0.1); 
    Td = 2* exp(-2*(ite/itemax)^2);  
    a=ceil((noThieves-1).*rand(noThieves,1))'; 
    
    for i=1:noThieves
         if (rand>=0.5)
                if rand>Pp 
                  xth(i,:)=gbest +(Td*(best(i,:)-xab(i,:))*rand+Td*(xab(i,:)-best(a(i),:))*rand)*sign(rand-0.50);
                else

                    for j=1:dim
                    xth(i,j)= Td*((ub(j)-lb(j))*rand+lb(j));
                    end

                end
         else 
                for j=1:dim
                    xth(i,j)=gbest(j) -(Td*(best(i,j)-xab(i,j))*rand+Td*(xab(i,j)-best(a(i),j))*rand)*sign(rand-0.50);          
                end
         end
    end
  
     for i=1:noThieves 
           fit(i,1)=fobj(xth(i,:));
           if and (~(xth(i,:)-lb<= 0),~(xth(i,:) - ub>= 0))            
               xab(i,:)=xth(i,:); 
                if fit(i)<fitness(i)
                     best(i,:) = xth(i,:); 
                     fitness(i)=fit(i); 
                end

                if fitness(i)<fit0
                   fit0=fitness(i);
                   gbest=best(i,:); 
                end
            end
     end

    ccurve(ite)=fit0; 
 
end

bestThieves=find(fitness==min(fitness)); 
gbestSol=best(bestThieves(1),:);  

fitness =fobj(gbestSol); 

end