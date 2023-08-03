clc
clear

str='E:\WMF\EXP\TSO8\';

for i=6    %K
    for g=1  %30
        for j=1  %image
            I1=imread([str,num2str(j),'.jpg']);
            I1=rgb2gray(I1);
            figure,imshow(I1);
            rgb = ind2rgb(gray2ind(I1,255),jet(255));
            figure,imshow(rgb,[],'border','tight')
                
        filepath=pwd;                                 
        cd('E:\WMF\1£¬K=2') 
                picname=['1,K=2-CS.jpg'];
        hold on 
        saveas(gcf,picname)
        cd(filepath) 
        end
    end                 
end