%% CTSO   A novel chaotic-driven Tuna Swarm Optimizer with Newton-Raphson method for parameter identification of three-diode equivalent circuit model of solar photovoltaic cells/ modules

function [Tuna1_fit,Tuna1,Convergence_curve]=CTSO(Particles_no,Max_iter,Low,Up,Dim,fobj)
Tuna1=zeros(1,Dim);   Tuna1_fit=inf;
pop=initialization(Particles_no,Dim,Up,Low);
Iter=0;
aa=0.7;  
z=0.05;
while Iter<Max_iter
    C=Iter/Max_iter;
    M=Tent1(Iter);
    qq=M(1,1);
    for i=1:size(pop,1)
        
        Flag4ub=pop(i,:)>Up;
        Flag4lb=pop(i,:)<Low;
        pop(i,:)=(pop(i,:).*(~(Flag4ub+Flag4lb)))+Up.*Flag4ub+Low.*Flag4lb;
        
        fitness(i)=fobj(pop(i,:));
        
        if fitness(i)<Tuna1_fit
            Tuna1_fit=fitness(i);  
            Tuna1=pop(i,:);
        end
    end
    if Iter==0
        fit_old=fitness;  C_old=pop;
    end
    
    for i=1:Particles_no
        if fit_old(i)<fitness(i)
            fitness(i)=fit_old(i); pop(i,:)=C_old(i,:);
        end
    end
    
    C_old=pop;  fit_old=fitness;
    t=(1-Iter/Max_iter)^(Iter/Max_iter);                   
    if rand<z
        pop(1,:)= (Up-Low)*rand+Low;
    else
        if  0.5>rand
            r1=rand;
            Beta=exp(r1*exp(3*cos(pi*((Max_iter-Iter+1)/Max_iter))))*(cos(2*pi*r1));
            if  C>rand
                pop(1,:)=qq.*(Tuna1+Beta*abs(Tuna1-pop(1,:)))+qq.*pop(1,:);
            else
                IndivRand=rand(1,Dim).*(Up-Low)+Low;
                pop(1,:)=qq.*(IndivRand+Beta*abs(IndivRand-pop(i,:)))+qq.*pop(1,:);
            end
        else
            TF = (rand>0.5)*2-1; 
            if 0.5<rand
                pop(1,:)=Tuna1+qq.*(Tuna1-pop(1,:))+TF.*t^2.*(Tuna1-pop(1,:));
            else
                pop(1,:) =TF.* t^2.*pop(1,:);
            end
        end
    end
    
    for i=2:Particles_no
        if rand<z    
            pop(i,:)= (Up-Low)*rand+Low;
        else
            if  0.5>rand
                r1=rand;
                Beta=exp(r1*exp(3*cos(pi*((Max_iter-Iter+1)/Max_iter))))*(cos(2*pi*r1));
                if  C>rand
                    pop(i,:)=qq.*(Tuna1+Beta*abs(Tuna1-pop(i,:)))+qq.*pop(i-1,:);
                else
                    IndivRand=rand(1,Dim).*(Up-Low)+Low;
                    pop(i,:)=qq.*(IndivRand+Beta*abs(IndivRand-pop(i,:)))+qq.*pop(i-1,:);
                end
            else
                TF = (rand>0.5)*2-1;
                if 0.5<rand
                    pop(i,:)=Tuna1+qq.*(Tuna1-pop(i,:))+TF*t^2.*(Tuna1-pop(i,:));
                else
                    pop(i,:) = TF*t^2.*pop(i,:);
                end
            end
        end
    end
    
    Iter=Iter+1;
    Convergence_curve(Iter)=Tuna1_fit;
    
end




