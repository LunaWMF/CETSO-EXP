close all
clear
clc

for number=1:30
    for i=1
        str='E:\WMF\Gray\';
        f0 = imread([str,num2str(i),'-Gray.jpg']);
        str='E:\WMF\\SymmetricCross-TSO8\';
        f1 = imread([str,num2str(i),'SymmetricCross-TSO8,K=7,',num2str(number),'.jpg']);  
        [PSNR, RMSE]   =    PSNR(f0, f1);   fprintf('%d',PSNR);     fprintf('\t'); 
        [mssim, ssim_map] = SSIM(f0, f1);   fprintf('%d',mssim);    fprintf('\t');
        [FSIM, FSIMc]  =    FSIM(f0, f1);   fprintf('%d',FSIM);     fprintf('\n');
        clear  
    end
end
