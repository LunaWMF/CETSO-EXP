%% grey entropy
function [g] = Grey2(I,thresh)
thresh = floor(thresh);
thresh = sort(thresh);
[r,c]=size(I);
 n=max(size(thresh));
 h = imhist(I);
 p= h./(r*c);
 B=size(1,n+1);
 HH=size(1,n+1);
 if(n==1)
     L1 = 1:thresh;
     L1 = L1';
     B0 = sum(L1.*h(L1));
     H= -1.*(h(L1).*L1./B0).*log(L1./B0);
     H(isnan(H))=0;
     HH(1)=sum(H);
     L2 = thresh+1:256;
     L2 = L2';
     B1 = sum(L2.*h(L2));
     H= -1.*(h(L2).*L2./B1).*log(L2./B1);
     H(isnan(H))=0;
     HH(2)=sum(H);     
 else
     for i=1:n+1
         if(i==1)
            L = 1:thresh(i);
            L = L';
            B(i) =sum(L.*h(L));
            H = -1.*(h(L).*L./B(i)).*log(L./B(i));
            H(isnan(H))=0;
            HH(i)=sum(H);
         elseif(i==n+1)
             L = thresh(i-1):256;
             L = L';
             B(i) =sum(L.*h(L));
             H = -1.*(h(L).*L./B(i)).*log(L./B(i));
             H(isnan(H))=0;
             HH(i)=sum(H);
         else
             L = thresh(i-1):thresh(i);
             L = L';
             B(i) =sum(L.*h(L));
             H = -1.*(h(L).*L./B(i)).*log(L./B(i));
             H(isnan(H))=0;
             HH(i)=sum(H);
         end
     end
 end
g=sum(HH);
end

