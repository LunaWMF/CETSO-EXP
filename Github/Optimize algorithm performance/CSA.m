function [fmin0,gPosition,cg_curve]=CSA(searchAgents,iteMax,lb,ub,dim,fobj)

if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
 
cg_curve=zeros(1,iteMax);

chameleonPositions=initialization(searchAgents,dim,ub,lb);

fit=zeros(searchAgents,1);

for i=1:searchAgents
     fit(i,1)=fobj(chameleonPositions(i,:));
end

fitness=fit; 
 
[fmin0,index]=min(fit);

chameleonBestPosition = chameleonPositions;
gPosition = chameleonPositions(index,:); 

v=0.1*chameleonBestPosition;
v0=0.0*v;

rho=1.0; 
gamma=2.0; 
alpha = 4.0;  
beta=3.0; 
 
for t=1:iteMax
    
    a = 2590*(1-exp(-log(t))); 
    omega=(1-(t/iteMax))^(rho*sqrt(t/iteMax)) ; 
    p1 = 2* exp(-2*(t/iteMax)^2);  
    p2 = 2/(1+exp((-t+iteMax/2)/100)) ;
    mu= gamma*exp(-(alpha*t/iteMax)^beta) ;
    ch=ceil(searchAgents*rand(1,searchAgents));

    for i=1:searchAgents  
        if rand>=0.1
           chameleonPositions(i,:)= chameleonPositions(i,:)+ p1*(chameleonBestPosition(ch(i),:)-chameleonPositions(i,:))*rand()+... 
           + p2*(gPosition -chameleonPositions(i,:))*rand();
        else 
        for j=1:dim
            chameleonPositions(i,j)=   gPosition(j)+mu*((ub(j)-lb(j))*rand+lb(j))*sign(rand-0.50) ;
        end 
        end   
    end       

     for i=1:searchAgents
        v(i,:)= omega*v(i,:)+ p1*(chameleonBestPosition(i,:)-chameleonPositions(i,:))*rand +.... 
               + p2*(gPosition-chameleonPositions(i,:))*rand;        

         chameleonPositions(i,:)=chameleonPositions(i,:)+(v(i,:).^2 - v0(i,:).^2)/(2*a);
     end
     v0=v;
  
     for i=1:searchAgents
         if chameleonPositions(i,:)<lb
            chameleonPositions(i,:)=lb;
         elseif chameleonPositions(i,:)>ub
            chameleonPositions(i,:)=ub;
         end
     end

    for i=1:searchAgents
        ub_=sign(chameleonPositions(i,:)-ub)>0;   
        lb_=sign(chameleonPositions(i,:)-lb)<0;
        chameleonPositions(i,:)=(chameleonPositions(i,:).*(~xor(lb_,ub_)))+ub.*ub_+lb.*lb_;  
        fit(i,1)=fobj (chameleonPositions(i,:)) ;
          if fit(i)<fitness(i)
             chameleonBestPosition(i,:) = chameleonPositions(i,:);
             fitness(i)=fit(i); 
          end
    end

    [fmin,index]=min(fitness); 
    if fmin < fmin0
        gPosition = chameleonBestPosition(index,:); 
        fmin0 = fmin;
    end
   cg_curve(t)=fmin0; 

end
 
ngPosition=find(fitness== min(fitness)); 
g_best=chameleonBestPosition(ngPosition(1),:); 
fmin0 =fobj(g_best);

end


%%
function pos=initialization(searchAgents,dim,u,l)
Boundary_no= size(u,2); 
if Boundary_no==1
    u_new=ones(1,dim)*u;
    l_new=ones(1,dim)*l;
else
     u_new=u;
     l_new=l;   
end
    for i=1:dim
        u_i=u_new(i);
        l_i=l_new(i);
        pos(:,i)=rand(searchAgents,1).*(u_i-l_i)+l_i;
    end
end
%%
function answer = get_orthonormal(m,n)
    if ( (nargin==2) && (m>n) && (isnumeric(m)*isnumeric(n)) )
    elseif ( nargin==1 && isnumeric(m) && length(m)==1 )
        n=m;
    else
       error('Incorrect Inputs. Please read help text in m-file.')
    end
    count=0;
    while (count==0)
        A=rand(m);
        B=A'*A ;
        [P,D] = eig(B);
        if ((P'*P - eye(m))>eps) 
            count=0;
        else
            answer=P(:,1:n) ;
            count=1;
        end
    end
end
%%
function [chameleonPositions]=rotation(chameleonPosition, searchAgents, dim)
for i=1:searchAgents      
          if (dim>2) 
              xmax=1;xmin=-1;
              th=round(xmin+rand(1,1)*(xmax-xmin));
              vec=get_orthonormal(dim,2);
              vecA=vec(:,1);
              vecB=vec(:,2);
              theta=(th*rand()*180)*(pi/180) ;
              Rot = RotMatrix(theta,vecA, vecB) ;
             if (theta~=0)
                V=[chameleonPosition(i,:) ]; 
                V_centre=mean(V,1); 
                Vc=V-ones(size(V,1),1)*V_centre;
                Vrc=[Rot*Vc']';
                Vr=Vrc+ones(size(V,1),1)*V_centre; 
                 chameleonPosition(i,:)=((Vr)/1); 
             end
         else
              xmax=1;xmin=-1;
              th=round(xmin+rand(1,1)*(xmax-xmin));
              theta=th*rand()*180*(pi/180);
              Rot = RotMatrix(theta);
              
               if (theta~=0)
                V=[chameleonPosition(i,:) ];  
                V_centre=mean(V,1); 
                Vc=V-ones(size(V,1),1)*V_centre; 

                Vrc=[Rot*Vc']'; 
                Vr=Vrc+ones(size(V,1),1)*V_centre; 
                chameleonPosition(i,:)=((Vr)/1);
               end
          end
end  
   chameleonPositions=chameleonPosition;
end
%%
function R = RotMatrix(alpha, u, v)

if numel(alpha) ~= 1
   error('JSimon:RotMatrrix:BadInput1', ...
      'Angle of rotation must be a scalar.');
end

s = sin(alpha);
c = cos(alpha);

switch nargin
   case 1
      R = [c, -s;  s, c];
   case 2
      if numel(u) ~= 3
         error('JSimon:RotMatrrix:BadAxis2D', ...
            '3D: Rotation axis must have 3 elements.');
      end
      u = u(:);
      u = u ./ sqrt(u.' * u);
      x  = u(1);
      y  = u(2);
      z  = u(3);
      mc = 1 - c;
      R  = [c + x * x * mc,      x * y * mc - z * s,   x * z * mc + y * s; ...
            x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s; ...
            x * z * mc - y * s,  y * z * mc + x * s,   c + z * z .* mc];
   case 3
      n = numel(u);
      if n ~= numel(v)
         error('JSimon:RotMatrrix:BadAxes3D', ...
            'ND: Axes to define plane of rotation must have the same size.');
      end
      u = u(:);
      u = u ./ sqrt(u.' * u);
      v = v(:);
      v = v - (u.' * v) * u;
      v = v ./ sqrt(v.' * v);
      R = eye(n) + ...
         (v * u.' - u * v.') * s + ...
         (u * u.' + v * v.') * (c - 1);
      
   otherwise
      error('JSimon:RotMatrrix:BadNInput', ...
            '1 to 3 inputs required.');
end

end

