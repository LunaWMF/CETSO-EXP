close all
clear
clc
str='E:\DATA\512-340\';
for i=1:6
    I=imread([str,num2str(i),'.jpg']);
    I=rgb2gray(I); 
    figure,imshow(I,[],'border','tight');
    colormap(gray);
      
        filepath=pwd;                                         
        cd('E:\DATA\512-340')
        picname=[num2str(i) ,'-Gray.jpg'];
        hold on 
        saveas(gcf,picname)
        cd(filepath)
        
end