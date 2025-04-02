function [wa]=solve_linear_equation(x,vx,vy,vz,WEIGHT)
% Solving a system of n * m (n<m) linear equations with special solutions containing 0.
epsilon=1e-2;
L_x=length(x);
wa=zeros(L_x,1); 
v1=zeros(L_x,L_x-1); v2=zeros(L_x,L_x-2);
v3=zeros(L_x,L_x-3); v4=zeros(L_x,L_x-4);
v5=zeros(L_x,L_x-5); v6=zeros(L_x,L_x-6);
v7=zeros(L_x,L_x-7); v8=zeros(L_x,L_x-8);
v0=eye(L_x);
 for k=1:L_x-1
     v1(:,k)=sum(v0(:,k)).*v0(:,k+1)-sum(v0(:,k+1)).*v0(:,k);
 end
 v1=v1./max(norm(v1,2),1e-20);
 for k=1:L_x-2
     sum1=sum(v1(:,k).*vx);
     sum2=sum(v1(:,k+1).*vx);
     if abs(sum1)<epsilon && abs(sum2)<epsilon 
         v2(:,k)=v1(:,k); 
     else
         v2(:,k)=sum1.*v1(:,k+1)-sum2.*v1(:,k); 
     end
 end
 v2=v2./max(norm(v2,2),1e-20);
 for k=1:L_x-3
     sum1=sum(v2(:,k).*vy);
     sum2=sum(v2(:,k+1).*vy);
     if abs(sum1)<epsilon && abs(sum2)<epsilon 
         v3(:,k)=v2(:,k);       
     else
         v3(:,k)=sum1.*v2(:,k+1)-sum2.*v2(:,k);  
     end
 end
 v3=v3./max(norm(v3,2),1e-20);
 for k=1:L_x-4
     sum1=sum(v3(:,k).*vz);
     sum2=sum(v3(:,k+1).*vz);
     if abs(sum1)<epsilon && abs(sum2)<epsilon 
         v4(:,k)=v3(:,k);         
     else 
         v4(:,k)=sum1.*v3(:,k+1)-sum2.*v3(:,k);
     end
 end
 v4=v4./max(norm(v4,2),1e-20);
 for k=1:L_x-5
     sum1=sum(v4(:,k).*vx.*vx);
     sum2=sum(v4(:,k+1).*vx.*vx);
     if abs(sum1)<epsilon && abs(sum2)<epsilon 
         v5(:,k)=v4(:,k); 
     else
         v5(:,k)=sum1.*v4(:,k+1)-sum2.*v4(:,k);
     end
 end
 v5=v5./max(norm(v5,2),1e-20);
 for k=1:L_x-6
     sum1=sum(v5(:,k).*vy.*vy);
     sum2=sum(v5(:,k+1).*vy.*vy);
     if abs(sum1)<epsilon && abs(sum2)<epsilon 
         v6(:,k)=v5(:,k);
     else
         v6(:,k)=sum1.*v5(:,k+1)-sum2.*v5(:,k);
     end
 end
 v6=v6./max(norm(v6,2),1e-20);
 for k=1:L_x-7
     sum1=sum(v6(:,k).*vz.*vz);
     sum2=sum(v6(:,k+1).*vz.*vz);
     if abs(sum1)<epsilon && abs(sum2)<epsilon 
         v7(:,k)=v6(:,k);
     else
         v7(:,k)=sum1.*v6(:,k+1)-sum2.*v6(:,k);
     end
 end
 v7=v7./max(norm(v7,2),1e-20);
 for k=1:L_x-8
     sum1=sum(v7(:,k).*x);
     sum2=sum(v7(:,k+1).*x);
     if abs(sum1)<epsilon && abs(sum2)<epsilon 
         v8(:,k)=v7(:,k);        
     else
         v8(:,k)=sum1.*v7(:,k+1)-sum2.*v7(:,k); 
     end
 end
 if max(abs(v8))>0
 v8(:,1)=v8(:,1)./max(norm(v8,2),1e-20);
 p_i=find(v8(:,1)<0);
 afa_p1=min(-WEIGHT(p_i)./v8(p_i,1));
 wa1=(WEIGHT+afa_p1.*v8(:,1));
 p_i=find(v8(:,1)>0);
 afa_p2=max(-WEIGHT(p_i)./v8(p_i,1));
 wa2=(WEIGHT+afa_p2.*v8(:,1));
 if rand<afa_p2/(afa_p2-afa_p1)
     wa=wa1;
 else
     wa=wa2;
 end
 else 
     wa=WEIGHT;
 end