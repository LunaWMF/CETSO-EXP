function [Iout]=YuZhiplot(I,thresh)
thresh = sort(thresh);
Iout=zeros(size(I));
n=max(size(thresh));
if(n==1)
Iout(I<thresh(1))=1;
Iout(I>=thresh(1))=2;
else
    for i=1:n+1
        if(i==1)
            Iout(I<thresh(i))=i;
         elseif(i==n+1)
             Iout(I>=thresh(i-1))=i;
        else
             Iout(I<thresh(i)&I>=thresh(i-1))=i;
         end
    end
end
imshow(Iout,[]);
end