close all
clear
clc
str='E:\DATA\512-340\';

    for m=1
        I=imread([str,num2str(m),'.jpg']); % 
        I=rgb2gray(I);         
        row=size(I,1);
        column=size(I,2);
        Z=zeros(1,256);
        for i=1:row
            for j=1:column
                k=I(i,j);
                Z(k+1)=Z(k+1)+1;
            end
        end
        figure;
        
        h=bar(Z);
        set(h,'FaceColor','[0.6,0.54,0.635]');
        axis tight;
        
    end
