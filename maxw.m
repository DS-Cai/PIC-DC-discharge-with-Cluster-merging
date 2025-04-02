function x=maxw(T_eV,N)
x=zeros(N,1);
% 通用气体常数 J/(mol K)
R_M=8.314;
% 玻尔兹曼常数
k_B=1.3807E-23;
% 电子电量 C
e=1.602176487E-19;
%电子质量 kg
m_e=9.1093822E-31;
T_K=e*T_eV/k_B;
v=1:2e6;
sumfv=zeros(2e6,1);
f_v=(m_e./2./pi()./k_B./T_K).^(3/2).*4.*pi().*v.*v.*exp(-m_e.*v.*v./2./k_B./T_K);
sumfv(1)=f_v(1);
for i=2:2e6
    sumfv(i)=sumfv(i-1)+f_v(i);
end
sumfv=sumfv./max(sumfv);
% U=rand(N,1);
U=1:N;
U=U./(N+1);
for i=1:N
x(i)=v(min(find(sumfv>U(i))));
end
% vx=zeros(20,1);
% nx=zeros(20,1);
% for i=1:20
%     vx(i)=(i-0.5)*1e5;
%     nx(i)=length(find(x>(i-1)*1e5 & x<(i)*1e5));
% end
% scatter(vx,nx)

