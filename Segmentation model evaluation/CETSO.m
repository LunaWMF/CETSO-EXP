
function [Tuna1_fit,Tuna1,Convergence_curve]=CETSO(Particles_no,Max_iter,Low,Up,Dim,fobj)
Tuna1=zeros(1,Dim);   
Tuna1_fit=inf;
pop=initialization(Particles_no,Dim,Up,Low);
Iter=0;
aa=0.7; 
z=0.05;
while Iter<Max_iter 
    C=Iter/Max_iter;
    a1=aa+(1-aa)*cos(C*pi/2);  
    a2=(1-aa)-(1-aa)*cos(C*pi/2);
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
    t=cos((Iter*pi)./Max_iter*25);
    if rand<z
        pop(1,:)= (Up-Low)*rand+Low;
    else
        if  0.5>rand
            b=rand; 
            Beta=exp(b*exp(3*cos(pi*((Max_iter-Iter+1)/Max_iter))))*(cos(2*pi*b));
            if  C>rand
                pop(1,:)=a1.*(Tuna1+Beta*abs(Tuna1-pop(1,:)))+a2.*pop(1,:);
            else
                IndivRand=rand(1,Dim).*(Up-Low)+Low;
                pop(1,:)=a1.*(IndivRand+Beta*abs(IndivRand-pop(i,:)))+a2.*pop(1,:);
            end
        else
            TF = (rand>0.5)*2-1; 
            if 0.5<rand
                pop(1,:)=Tuna1+rand(1,Dim).*(Tuna1-pop(1,:))+TF.*t^2.*(Tuna1-pop(1,:));
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
                b=rand;
                Beta=exp(b*exp(3*cos(pi*((Max_iter-Iter+1)/Max_iter))))*(cos(2*pi*b));
                if  C>rand
                    pop(i,:)=a1.*(Tuna1+Beta*abs(Tuna1-pop(i,:)))+a2.*pop(i-1,:);
                else
                    
                    IndivRand=rand(1,Dim).*(Up-Low)+Low;
                    pop(i,:)=a1.*(IndivRand+Beta*abs(IndivRand-pop(i,:)))+a2.*pop(i-1,:);
                end
            else
                TF = (rand>0.5)*2-1;
                if 0.5<rand
                    pop(i,:)=Tuna1+rand(1,Dim).*(Tuna1-pop(i,:))+TF*t^2.*(Tuna1-pop(i,:));
                else
                    pop(i,:) = TF*t^2.*pop(i,:);
                end
            end
        end
    end
    
    b = 1-C.^2;
    r2 = C.^2;
    cauchy = rand();  
    gauss = rand();
    pop(i,:) = pop(1,:).*(1+b.*cauchy+r2.*gauss);
    
    Iter=Iter+1;
    Convergence_curve(Iter)=Tuna1_fit;
end




