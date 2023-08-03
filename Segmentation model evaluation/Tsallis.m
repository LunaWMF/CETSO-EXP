%% Tsallis
function [J]=Tsallis(I,thresh)
[r,c]=size(I);
t= round(thresh);
n=max(size(thresh));
Ihist=imhist(I)/(r*c);
Ihistcum=cumsum(Ihist);
q=0.8;

P=zeros(1,n+1);
u=zeros(1,n+1);
deta=zeros(1,n+1);
pp=zeros(n+1,256);
L=1:256;
R=zeros(1,256);
if(n==1)
    P(1)=Ihistcum(t);
    P(2)=Ihistcum(end)-Ihistcum(t);
    u(1)=sum(Ihist(1:t).*(1:t)')/P(1);
    u(2)=sum(Ihist(t+1:256).*(t+1:256)')/P(2);
    u(isnan(u))=0;
    deta(1)=sum(((1:t)-u(1).*ones(size(1:t))).^2)/P(1);
    deta(2)=sum( ((t+1:256)-u(2).*ones(size(t+1:256))).^2)/P(2);
    deta(isnan(deta))=0;
    pp(1,:)=normpdf(L,u(1),deta(1)^0.5);
    pp(2,:)=normpdf(L,u(2),deta(2)^0.5);
    R=P(1).*pp(1,:)+P(2).*pp(2,:);
else
    for i=1:n+1
        if(i==1)
            P(i)=Ihistcum(t(i));
            u(i)=sum(Ihist(1:t(i)).*(1:t(i))')/sum(Ihist(1:t(i)));
            u(isnan(u))=0;
            deta(i)=sum(((1:t(i))-u(i).*ones(size(1:t(i)))).^2)/sum(Ihist(1:t(i)));
            deta(isnan(deta))=0;
            pp(i,:)=normpdf(L,u(i),deta(i)^0.5);
        elseif(i==n+1)
          P(i)=Ihistcum(end)-Ihistcum(t(i-1)+1);
          u(i)=sum(Ihist(t(i-1)+1:256).*(t(i-1)+1:256)')/sum(Ihist(t(i-1)+1:256));
          u(isnan(u))=0;
          deta(i)=sum( ((t(i-1)+1:256)-u(i).*ones(size(t(i-1)+1:256))).^2)/sum(Ihist(t(i-1)+1:256));
          deta(isnan(deta))=0;
          pp(i,:)=normpdf(L,u(i),deta(i)^0.5);
        else
            P(i)=Ihistcum(t(i))-Ihistcum(t(i-1)+1);
            u(i)=sum(Ihist(t(i-1)+1:t(i)).*(t(i-1)+1:t(i))')/sum(Ihist(t(i-1)+1:t(i)));
            u(isnan(u))=0;
            deta(i)=sum(((t(i-1)+1:t(i))-u(i).*ones(size(t(i-1)+1:t(i)))).^2)/sum(Ihist(t(i-1)+1:t(i)));
            deta(isnan(deta))=0;
            pp(i,:)=normpdf(L,u(i),deta(i)^0.5);
        end
        R=R+P(i).*pp(i,:);
    end
end
temp1=Ihist'.*(R./Ihist').^q;
temp1(isnan(temp1))=0;
temp2=R.*(Ihist'./R).^q;
temp2(isnan(temp2))=0;
temp=temp1+temp2;
J=(sum(temp)-2)/(q-1);    
end