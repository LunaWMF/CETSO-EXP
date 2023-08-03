%%White Shark Optimizer
function [fmin0,gbest,ccurve]=WSO(whiteSharks,itemax,lb,ub,dim,fobj)
ccurve=zeros(1,itemax);
WSO_Positions=initialization(whiteSharks,dim,ub,lb);
v=0.0*WSO_Positions; 
fit=zeros(whiteSharks,1);

for i=1:whiteSharks
     fit(i,1)=fobj(WSO_Positions(i,:));
end
fitness=fit; 
[fmin0,index]=min(fit);

wbest = WSO_Positions; 
gbest = WSO_Positions(index,:); 

    fmax=0.75; 
    fmin=0.07; 
    tau=4.11;  
       
    mu=2/abs(2-tau-sqrt(tau^2-4*tau));

    pmin=0.5;
    pmax=1.5;
    a0=6.250;  
    a1=100;
    a2=0.0005;
    
for ite=1:itemax
    mv=1/(a0+exp((itemax/2.0-ite)/a1)); 
    s_s=abs((1-exp(-a2*ite/itemax))) ;
    p1=pmax+(pmax-pmin)*exp(-(4*ite/itemax)^2);
    p2=pmin+(pmax-pmin)*exp(-(4*ite/itemax)^2);
     nu=floor((whiteSharks).*rand(1,whiteSharks))+1;

     for i=1:size(WSO_Positions,1)
           rmin=1; rmax=3.0;
          rr=rmin+rand()*(rmax-rmin);
          wr=abs(((2*rand()) - (1*rand()+rand()))/rr);       
          v(i,:)=  mu*v(i,:) +  wr *(wbest(nu(i),:)-WSO_Positions(i,:));      
     end
 
     for i=1:size(WSO_Positions,1)
        f =fmin+(fmax-fmin)/(fmax+fmin);
        a=sign(WSO_Positions(i,:)-ub)>0;
        b=sign(WSO_Positions(i,:)-lb)<0;
        wo=xor(a,b);
            if rand<mv
                WSO_Positions(i,:)=  WSO_Positions(i,:).*(~wo) + (ub.*a+lb.*b);
            else   
                WSO_Positions(i,:) = WSO_Positions(i,:)+ v(i,:)/f;  
            end
    end 
for i=1:size(WSO_Positions,1)
        for j=1:size(WSO_Positions,2)
            if rand<s_s      
                
             Dist=abs(rand*(gbest(j)-1*WSO_Positions(i,j)));
             
                if(i==1)
                    WSO_Positions(i,j)=gbest(j)+rand*Dist*sign(rand-0.5);
                else    
                    WSO_Pos(i,j)= gbest(j)+rand*Dist*sign(rand-0.5);
                    WSO_Positions(i,j)=(WSO_Pos(i,j)+WSO_Positions(i-1,j))/2*rand;
                end   
            end
         
        end       
end
 
    for i=1:whiteSharks 
           if WSO_Positions(i,:)>=lb & WSO_Positions(i,:)<=ub
              fit(i)=fobj(WSO_Positions(i,:));   
            if fit(i)<fitness(i)
                 wbest(i,:) = WSO_Positions(i,:); 
                 fitness(i)=fit(i);  
            end
            if (fitness(i)<fmin0)
               fmin0=fitness(i);
               gbest = wbest(index,:);
            end 
            
        end
    end
  
end 
end

function pos=initialization(whiteSharks,dim,ub_,lb_)
BoundNo= size(ub_,1); 
if BoundNo==1
    pos=rand(whiteSharks,dim).*(ub_-lb_)+lb_;
end
if BoundNo>1
    for i=1:dim
        ubi=ub_(i);
        lbi=lb_(i);
        pos(:,i)=rand(whiteSharks,1).*(ubi-lbi)+lbi;
    end
end
end