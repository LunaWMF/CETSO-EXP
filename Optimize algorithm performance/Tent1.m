function [x] = Tent1(maxIter)
 a=0.7;
 x(1)=rand;
     for i=1:maxIter-1
         if x(i)<a
             x(i+1)=x(i)/a;
         end
         if x(i)>=a
             x(i+1)=(1-x(i))/a;
         end
     end
 end


 