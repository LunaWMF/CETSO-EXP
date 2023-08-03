

close all
clear
clc
rng('default');
str='E:\WMF\Code\Segmentation model evaluation\DATA_image\';

Particles_no=30; 
Max_iter=100;

for number=1:30  %   :30
    for i=3  %K      :2:7 
        for j=1  %IMAGE  :6
            I=imread([str,num2str(j),'.jpg']);
            I=rgb2gray(I);        
            dim = i;
            lb = ones(1,dim);
            ub = 255.*ones(1,dim);
        
            fobj =@(thresh)EXP(I,thresh);% EXP  Grey  Tsallis  MinCross   SymCross
           
%             [Best_score1,Best_pos1,TSO_cg_curve1]=TSO(Particles_no,Max_iter,lb,ub,dim,fobj);
%             Best_pos1 = round(Best_pos1);
%             Best_Thresh1 = Best_pos1;
%             Best_Thresh1 = sort(Best_Thresh1);
%             [Iout1]=YuZhiplot(I,Best_Thresh1);
%             figure,imshow(Iout1,[],'border','tight');
%              display(['TSO,',num2str(j) ,',',num2str(i),',',num2str(number) ,',', num2str(Best_pos1),',', num2str(-Best_score1)]);
             
             [Best_score2,Best_pos2,CETSO_cg_curve2]=CETSO(Particles_no,Max_iter,lb,ub,dim,fobj);
            Best_pos2 = round(Best_pos2);
            Best_Thresh2 = Best_pos2;
            Best_Thresh2 = sort(Best_Thresh2);
            [Iout2]=YuZhiplot(I,Best_Thresh2);
            figure,imshow(Iout2,[],'border','tight');
             display(['CETSO,',num2str(j) ,',',num2str(i),',',num2str(number),',' , num2str(Best_pos2),',',num2str(-Best_score2)]);
             
            filepath=pwd;                                         
            cd('E:\WMF\Code\Segmentation model evaluation\EXP_Seg') 
             picname=[num2str(j),'EXP-CETSO,K=',num2str(i) ,',',num2str(number),'.jpg'];
            hold on 
            saveas(gcf,picname)
            cd(filepath) 

        end
%           fprintf('\n');
    end
end

