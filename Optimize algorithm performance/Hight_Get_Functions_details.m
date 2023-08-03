


function [lb,ub,dim,fobj] = Hight_Get_Functions_details(F)
    d=30;   
    switch F
        case 'F1'
            fobj = @F1;lb=-100;ub=100;dim=d;
        case 'F2'
            fobj = @F2;lb=-1;ub=1;dim=d;
        case 'F3'
            fobj = @F3;lb=-100;ub=100;dim=d;
        case 'F4'
            fobj = @F4;lb=-10;ub=10;dim=d;
        case 'F5'
            fobj = @F5;lb=-100;ub=100;dim=d;
        case 'F6'
            fobj = @F6;lb=-50;ub=50;dim=d;
        case 'F7'
            fobj = @F7;lb=-10;ub=10;dim=d;
        case 'F8'
            fobj = @F8;lb=-10;ub=10;dim=d;
    end
end

function o = F1(x)
    dim=size(x,2);
    o=0;
    for i=1:dim
        o=o+sum(x(1:i))^2;
    end
end

function o = F2(x)
    dim=size(x,2);
    o=0;
    for i=1:dim
         o=o+(abs(x(i))^(i+1));
    end
end

function o = F3(x)
    dim=size(x,2);
    o=0;
    for i=2:dim
    o=sum(x(i))^2;
    end
    o=x(1)^2+10^6*o;
end

function o = F4(x)
    dim=size(x,2);
    o=0;
    for i=1:dim
        o=o+i*x(i)^2;
    end
end

function o = F5(x)
    o=1-cos(2*pi*sqrt(sum(x.^2)))+0.1*sqrt(sum(x.^2));
end

function o = F6(x)
    dim=size(x,2);
    o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

function o = F7(x)
    dim=size(x,2);
    o=sum((x(1:dim-1)-1).^2.*(1+(sin(3*pi.*(x(2:dim))))))+sin(3*pi*x(1))^2+...
    abs((x(dim)-1)).*(1+(((sin(2*pi.*(x(dim))))).^2));
end

function o = F8(x)
    o=sum(abs(x.*sin(x)+0.1.*x));
end

