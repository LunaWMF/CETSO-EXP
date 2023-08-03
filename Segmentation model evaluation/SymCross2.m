%% symmetric cross Entropy
function [g] = SymCross2(I,thresh)
thresh = floor(thresh);
thresh = sort(thresh);
[r,c]=size(I);
 n=max(size(thresh));
 h = imhist(I);
 p= h./(r*c);
 u=size(1,n+1);
 HH=size(1,n+1);
 if(n==1)
     L1 = 1:thresh;
     L1 = L1';
     u0 = sum(L1.*p(L1))/sum(p(L1));
     H= h(L1).*(L1.*log(L1./u0)+u0.*log(u0./L1));
     H(isnan(H))=0;
     HH(1)=sum(H);
     L2 = thresh+1:256;
     L2 = L2';
     u1 = sum(L2.*p(L2))/sum(p(L2));
     H=h(L2).*(L2.*log(L2./u1)+u1.*log(u1./L2));
     H(isnan(H))=0;
     HH(2)=sum(H);     
 else
     for i=1:n+1
         if(i==1)
            L = 1:thresh(i);
            L = L';
            u(i) = sum(L.*p(L))/sum(p(L));
            H = h(L).*(L.*log(L./u(i))+u(i).*log(u(i)./L));
            H(isnan(H))=0;
            HH(i)=sum(H);
         elseif(i==n+1)
             L = thresh(i-1):256;
             L = L';
             u(i) = sum(L.*p(L))/sum(p(L));
             H = h(L).*(L.*log(L./u(i))+u(i).*log(u(i)./L));
             H(isnan(H))=0;
             HH(i)=sum(H);
         else
             L = thresh(i-1):thresh(i);
             L = L';
             u(i) = sum(L.*p(L))/sum(p(L));
             H = h(L).*(L.*log(L./u(i))+u(i).*log(u(i)./L));
             H(isnan(H))=0;
             HH(i)=sum(H);
         end
     end
 end
g=sum(HH);
end

