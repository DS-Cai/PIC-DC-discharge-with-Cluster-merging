clear all; clc;
% T_eV_=1;
% T_K=e*T_eV/k_B;
% syms v;
% f_v=v*(m_e/2/pi()/k_B/T_K)^(3/2)*4*pi()*v*v*exp(-m_e*v*v/2/k_B/T_K);
%通用气体常数 J/(mol K)
R_M=8.314;
%玻尔兹曼常数
k_B=1.3807E-23;
%电子电量 C
e=1.602176487E-19;
%电子质量 kg
m_e=9.1093822E-31;
NN=2e4; 
x=rand(NN,1);
WEIGHT=100.*ones(NN,1);
DX=0.01;
N_G=1/DX+1;
test_minsize=[8 12 16 20];
for k=1:4
min_size=test_minsize(k);
vx1=maxw(1,NN/2);
vx2=-vx1;
vx=[vx1;vx2];
vy=vx;
vz=vx;
[x2,vx2,vy2,vz2,WEIGHT2]=merging1(x,vx,vy,vz,WEIGHT,DX,N_G,min_size);
DN=50;
nvx=zeros(DN,2);
v_ave=zeros(DN,2); 
v_sqr=vx;
minvx=min(v_sqr);
dv=(max(v_sqr)-minvx)./DN;
for i=1:DN
    v_ave(i,1)=dv*(i-0.5)+minvx;
    ii=find(v_sqr>minvx+(i-1)*dv & v_sqr<=minvx+i*dv );
    if length(ii)
    nvx(i,1)=sum(WEIGHT(ii));
    else
        nvx(i,1)=0;
    end
end
v_sqr2=vx2;
minvx2=min(v_sqr2);
dv2=(max(v_sqr2)-minvx2)./DN;
for i=1:DN
    v_ave(i,2)=dv*(i-0.5)+minvx2;
    ii=find(v_sqr2>minvx2+(i-1)*dv2 & v_sqr2<=minvx2+i*dv2 );
    if length(ii)
    nvx(i,2)=sum(WEIGHT2(ii));
    else
        nvx(i,2)=0;
    end
end
figure
hold on
scatter(v_ave(:,1),nvx(:,1))
scatter(v_ave(:,2),nvx(:,2))
ave=[v_ave(:,1),nvx(:,1),v_ave(:,2),nvx(:,2)];
xlabel('Velocity'); % 设置横坐标标签
ylabel('Particle Number'); % 设置纵坐标标签
title(['Merge ' num2str(k)]); % 设置图表标题
grid on; % 显示网格
hold off
end