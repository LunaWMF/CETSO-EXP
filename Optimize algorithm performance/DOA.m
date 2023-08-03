% Dingo Optimization Algorithm (DOA) 
function [vMin,theBestVct,Convergence_curve]=DOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
P= 0.5;  
Q= 0.7; 
beta1= -2 + 4* rand(); 
beta2= -1 + 2* rand();  
naIni= 2;
naEnd= SearchAgents_no /naIni; 
na= round(naIni + (naEnd-naIni) * rand()); 
Positions=initialization(SearchAgents_no,dim,ub,lb);
 for i=1:size(Positions,1)
      Fitness(i)=fobj(Positions(i,:));      
 end
[vMin minIdx]= min(Fitness); 
theBestVct= Positions(minIdx,:); 
[vMax maxIdx]= max(Fitness); 
Convergence_curve=zeros(1,Max_iter);
Convergence_curve(1)= vMin;
survival= survival_rate(Fitness,vMin,vMax);  
t=0;

for t=1:Max_iter       
   for r=1:SearchAgents_no
      sumatory=0;
    if rand() < P  
           sumatory= Attack( SearchAgents_no, na, Positions, r );    
           if rand() < Q            
                 v(r,:)=  beta1 * sumatory-theBestVct; 
           else  
               r1= round(1+ (SearchAgents_no-1)* rand()); 
               v(r,:)= theBestVct + beta1*(exp(beta2))*((Positions(r1,:)-Positions(r,:))); 
           end  
    else 
        r1= round(1+ (SearchAgents_no-1)* rand());
        v(r,:)=   (exp(beta2)* Positions(r1,:)-((-1)^getBinary)*Positions(r,:))/2;
    end
    if survival(r) <= 0.3
         band=1; 
         while band 
           r1= round(1+ (SearchAgents_no-1)* rand());
           r2= round(1+ (SearchAgents_no-1)* rand());
           if r1 ~= r2 
               band=0;
           end
         end
              v(r,:)=   theBestVct + (Positions(r1,:)-((-1)^getBinary)*Positions(r2,:))/2;  
    end 
   
        Flag4ub=v(r,:)>ub;
        Flag4lb=v(r,:)<lb;
        v(r,:)=(v(r,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

     Fnew= fobj(v(r,:));

     if Fnew <= Fitness(r)
        Positions(r,:)= v(r,:);
        Fitness(r)= Fnew;
     end
     if Fnew <= vMin
         theBestVct= v(r,:);
         vMin= Fnew;
     end 
   end
   Convergence_curve(t+1)= vMin; 
   [vMax , ~]= max(Fitness);
    survival= survival_rate( Fitness, vMin, vMax);
end
end
 

function [ vAttack ] = vectorAttack( SearchAgents_no,na )
c=1; 
vAttack=[];
 while(c<=na)
    idx =round( 1+ (SearchAgents_no-1) * rand());
    if ~findrep(idx, vAttack)
        vAttack(c) = idx;
        c=c+1;
    end
 end
end

function [ o ] =  survival_rate(  fit, min, max )
    for i=1:size(fit,2)
         o(i)= (max-fit(i))/(max-min);
    end
end

function [ val] = getBinary( )
if rand() < 0.5
     val= 0;
else
     val=1;
end
end

function [ band ] = findrep( val, vector )
band= 0;
for i=1:size(vector, 2)
    if val== vector(i)
        band=1;
        break;
    end
end
end

function [ sumatory ] = Attack( SearchAgents_no, na, Positions, r )
sumatory=0;
vAttack= vectorAttack( SearchAgents_no, na );
     for j=1:size(vAttack,2)
           sumatory= sumatory + Positions(vAttack(j),:)- Positions(r,:);                       
     end 
     sumatory=sumatory/na;
end
