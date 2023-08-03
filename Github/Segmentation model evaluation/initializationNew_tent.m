function Positions=initializationNew_tent(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); 

    if Boundary_no==1
        for i = 1:SearchAgents_no
            Temp = Tent(dim).*(ub-lb)+lb;
            Temp(Temp>ub) = ub;
            Temp(Temp<lb) = lb;
            Positions(i,:)=Temp;
        end
    end

    if Boundary_no>1
       for j = 1:SearchAgents_no
           TentValue =  Tent(dim);
            for i=1:dim
                ub_i=ub(i);
                lb_i=lb(i);
                Positions(j,i)=TentValue(i).*(ub_i-lb_i)+lb_i;
                if(Positions(j,i)>ub_i)
                    Positions(j,i) = ub_i;
                end
                if(Positions(j,i)<lb_i)
                    Positions(j,i) = lb_i;
                end
            end
       end
    end
    
    
end




