function func_plot_con(func_name)

[lb,ub,dim,fobj] = Hight_Get_Functions_details(func_name);

switch func_name 
    case 'F1' 
        x=-100:2:100; y=x;      
    case 'F2' 
        x=-1:0.01:1;y=x;   
    case 'F3' 
        x=-100:2:100; y=x;     
    case 'F4' 
        x=-10:0.5:10; y=x; 
    case 'F5' 
        x=-100:2:100; y=x;
    case 'F6' 
        x=-5.12:0.1:5.12; y=x;
    case 'F7' 
        x=-10:0.5:10;  y=x;  
    case 'F8' 
        x=-10:0.1:10;   y=x;   
end    

L=length(x);
f=[];

for i=1:L
    for j=1:L
            f(i,j)=fobj([x(i),y(j)]);         
    end
end

surfc(x,y,f,'LineStyle','none');
colormap parula

end
