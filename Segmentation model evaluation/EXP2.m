%% exponential entropy
function [g] = EXP2(I,thresh)
thresh = floor(thresh);
thresh = sort(thresh);
[r,c]=size(I);
 n=max(size(thresh));
 h = imhist(I);
 p= h./(r*c);
 HH=zeros(1,n+1);
 w = zeros(1,n+1);
 if(n==1)
     L1 = 1:thresh;
     L1 = L1';
     w0 = sum(p(L1));
     H= p(L1).* exp(1 - p(L1)./w0)/w0;
     H(isnan(H))=0;
     HH(1)=sum(H);
     L2 = thresh+1:256;
     L2 = L2';
     w1 = sum(p(L2));
     H= p(L2).* exp(1 - p(L2)./w1)/w1;
     H(isnan(H))=0;
     HH(2)=sum(H);     
 else
     for i=1:n+1
         if(i==1)
            L = 1:thresh(i);
            L = L';
            w(i) = sum(p(L));
            H= p(L).* exp(1 - p(L)./w(i))/w(i);
            H(isnan(H))=0;
            HH(i)=sum(H);
         elseif(i==n+1)
             L = thresh(i-1):256;
             L = L';
            w(i) = sum(p(L));
            H= p(L).* exp(1 - p(L)./w(i))/w(i);
             H(isnan(H))=0;
             HH(i)=sum(H);
         else
             L = thresh(i-1):thresh(i);
             L = L';
             w(i) = sum(p(L));
             H= p(L).* exp(1 - p(L)./w(i))/w(i);
             H(isnan(H))=0;
             HH(i)=sum(H);
         end
     end
 end
g=sum(HH);
end

