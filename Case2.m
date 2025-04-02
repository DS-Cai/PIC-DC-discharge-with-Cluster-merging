close all;clear all;clc;
k=0.7;
wr=1.67387; wi=-0.392401;

% parameters
n_weit=10;
L=2*pi/k; dt=.02; nt=3000; ng=100; np=30000; weit=n_weit.*ones(np,1);
vb=1.0; xp1=1.0e-2; vp1=0.0; 
vt=0.3; % note: the normalization sqrt(2) will be found in randn()
wp=1; qm=-1;
q=wp^2/(qm*n_weit*np/L); rho_back=-q*n_weit*np/L; dx=L/ng;
EEk=zeros(4,nt); EEf=zeros(4,nt); lndE=zeros(4,nt);
% initial loading for the 2 Stream instability
xpp=linspace(0,L,np)';
vpp=vt*randn(np,1)+(1-2*mod([1:np]',2)).*vb;
vp=vpp; xp=xpp;
% vp=vt*randn(np,1); % randn is {exp[-(x-mu)^2/(2*sgm^2)]}/[sgm*sqrt(2*pi)]

% Perturbation
vp=vp+vp1*cos(k*xp);
xp=xp+xp1*cos(k*xp);
p=1:np;p=[p p];
ri=find(vp>0); bi=find(vp<=0);
% Main computational cycle
h = figure('Unit','Normalized','position',...
    [0.02 0.3 0.6 0.6],'DefaultAxesFontSize',15);
for it=1:nt
    % apply periodic bc on the particle positions
    xp=xp./L+10.0; xp=L.*(xp-floor(xp));
    
    % diagnosing
    if(mod(it,nt/4)==1)
        xr=[xp(ri)./L,vp(ri)./max(abs(vpp))];
        xb=[xp(bi)./L,vp(bi)./max(abs(vpp))];
        subplot(2,2,floor(4*it/nt)+1); 
        plot(xp(ri),vp(ri),'r.',xp(bi),...
            vp(bi),'b.','Markersize',2); 
        axis([0,L,-3*(abs(vt)+abs(vb)),3*(abs(vt)+abs(vb))]);
        title(['Without Merging, t=',num2str((it-1)*dt)]);
        xlabel('xp');ylabel('vp'); pause(0.2);
%         print(gcf, '-dpng', ['vp-x,t=',num2str(it*dt),'.png']);
    end
    
    % update xp
    xp=xp+vp*dt;
    
    % projection p->g 
    g1=floor(xp/dx-.5)+1;g=[g1;g1+1];
    fraz1=weit.*(1-abs(xp/dx-g1+.5));
    fraz=[fraz1;weit-fraz1];
    
    % apply bc on the projection
    out=(g<1);g(out)=g(out)+ng;
    out=(g>ng);g(out)=g(out)-ng;
    mat=sparse(p,g,fraz,np,ng);
    rho=full((q/dx)*sum(mat))'+rho_back;
    
    % computing fields, dE/dx
    Eg=zeros(ng,1);
    for j=1:ng-1
        Eg(j+1)=Eg(j)+(rho(j)+rho(j+1))*dx/2;
    end
    Eg(1)=Eg(ng)+rho(ng)*dx;
    Eg=Eg-mean(Eg);
    
    % projection q->p and update of vp
    vp=vp+mat*qm*Eg*dt./weit;
    
    EEk(it,1)=0.5*abs(q)*sum(weit.*vp.^2); % kinetic energy  
    EEf(it,1)=0.5*sum(Eg.^2)*dx; % potential energy
    t(it)=it*dt;
end
test_minsize=[8 12 16];
for k=1:3
% parameters
n_weit=10;
np=30000; weit=n_weit.*ones(np,1);
wp=1; qm=-1;  
q=wp^2/(qm*n_weit*np/L); rho_back=-q*n_weit*np/L; dx=L/ng;
% vpp=vt*randn(np,1)+(1-2*mod([1:np]',2)).*vb;
vp=vpp; xp=xpp;
% vp=vt*randn(np,1); % randn is {exp[-(x-mu)^2/(2*sgm^2)]}/[sgm*sqrt(2*pi)]
% Perturbation
vp=vp+vp1*cos(k*xp);
xp=xp+xp1*cos(k*xp);
constraint_num=test_minsize(k);
vy=vp; vz=vy;
[xp,vp,vy,vz,weit]=merging1(xp,vp,vy,vz,weit,dx,ng,constraint_num);
np=length(xp);p=1:np;p=[p p];
ri=find(vp>0); bi=find(vp<=0);
% Main computational cycle
h = figure('Unit','Normalized','position',...
    [0.02 0.3 0.6 0.6],'DefaultAxesFontSize',15);
title(['Merge ' num2str(k)]); % 设置图表标题
for it=1:nt
    % apply periodic bc on the particle positions
    xp=xp./L+10.0; xp=L.*(xp-floor(xp));
    
    % diagnosing
    if(mod(it,nt/4)==1)
        subplot(2,2,floor(4*it/nt)+1); 
        plot(xp(ri),vp(ri),'r.',xp(bi),...
            vp(bi),'b.','Markersize',2); 
        axis([0,L,-3*(abs(vt)+abs(vb)),3*(abs(vt)+abs(vb))]);
        title(['Merge ', num2str(k+1),', t=',num2str((it-1)*dt)]);
        xlabel('xp');ylabel('vp'); pause(0.2);
%         print(gcf, '-dpng', ['vp-x,t=',num2str(it*dt),'.png']);
    end
    
    % update xp
    xp=xp+vp*dt;
    
    % projection p->g 
    g1=floor(xp/dx-.5)+1;g=[g1;g1+1];
    fraz1=weit.*(1-abs(xp/dx-g1+.5));
    fraz=[fraz1;weit-fraz1];
    
    % apply bc on the projection
    out=(g<1);g(out)=g(out)+ng;
    out=(g>ng);g(out)=g(out)-ng;
    mat=sparse(p,g,fraz,np,ng);
    rho=full((q/dx)*sum(mat))'+rho_back;
    
    % computing fields, dE/dx
    Eg=zeros(ng,1);
    for j=1:ng-1
        Eg(j+1)=Eg(j)+(rho(j)+rho(j+1))*dx/2;
    end
    Eg(1)=Eg(ng)+rho(ng)*dx;
    Eg=Eg-mean(Eg);
    
    % projection q->p and update of vp
    vp=vp+mat*qm*Eg*dt./weit;
    
    EEk(it,k+1)=0.5*abs(q)*sum(weit.*vp.^2); % kinetic energy  
    EEf(it,k+1)=0.5*sum(Eg.^2)*dx; % potential energy
    t(it)=it*dt;
end
end
%%
figure
hold on
plot(t, EEk(:,1) + EEf(:,1), 'r', 'LineWidth',2);
plot(t, EEk(:,1), 'k', t, EEf(:,1), 'k:', 'LineWidth',2);
plot(t, EEk(:,2), 'm', t, EEf(:,2), 'm:', 'LineWidth',2);
plot(t, EEk(:,3), 'b', t, EEf(:,3), 'b:', 'LineWidth',2);
plot(t, EEk(:,4), 'g', t, EEf(:,4), 'g:', 'LineWidth',2);
hold off
% title(['(a) k=',num2str(k),', \omega_{theory}=',...
%     num2str(wr+1i*wi)],'fontsize',15);
xlabel('t'); ylabel('Energy');legend('E_k','E_e','E_{tot}');
