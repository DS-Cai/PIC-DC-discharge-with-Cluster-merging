% function [N_e,N_i,x_e,vx_e,vy_e,vz_e,x_i,vx_i,vy_i,vz_i]=collision_electron(xe,...
%      vxe,vye,vze,energy_index,N_e,N_i,x_e,vx_e,vy_e,vz_e,x_i,vx_i,vy_i,vz_i,E_MASS,AR_MASS,E_EXC_TH,EV_TO_J,E_ION_TH,sigma_e,RMB)
%  [N_e,N_i,x_e,vx_e,vy_e,vz_e,x_i,vx_i,vy_i,vz_i]= collision_electron(x_e(k+1), vx_e(k+1),...
%             vy_e(k+1), vz_e(k+1),energy_index,N_e,N_i,x_e,vx_e,vy_e,vz_e,x_i,vx_i,vy_i,vz_i,E_MASS,AR_MASS,E_EXC_TH,EV_TO_J,E_ION_TH,sigma_e,RMB);
%         N_e_coll = N_e_coll + 1;
xe=x_e(k+1); vxe=vx_e(k+1); vye=vy_e(k+1); vze=vz_e(k+1); 
%% 变量定义
F1 = E_MASS / (E_MASS + AR_MASS); F2 = AR_MASS / (E_MASS + AR_MASS);
t0 = 0; t1 = 0; t2 = 0; rnd = 0; g = 0; g2 = 0; gx = 0; gy = 0; gz = 0;
wx = 0; wy = 0; wz = 0; theta = 0; phi = 0; chi = 0; eta = 0; chi2 = 0;
eta2 = 0; sc = 0; cc = 0; se = 0; ce = 0; st = 0; ct = 0; sp = 0; cp = 0;
energy = 0; e_sc = 0; e_ej = 0;
%% calculate relative velocity before collision & velocity of the centre of mass
gx = vxe; gy = vye; gz = vze;
g = sqrt(gx * gx + gy * gy + gz * gz);
wx = F1 * vxe; wy = F1 * vye; wz = F1 * vze;
%% find Euler angles
if gx == 0
    theta = 0.5 * pi;
else
    theta = atan2(sqrt(gy * gy + gz * gz), gx);
end
if gy == 0
    if gz > 0
        phi = 0.5 * pi;
    else
        phi = -0.5 * pi;
    end 
else
    phi = atan2(gz, gy);
end
st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);
%% choose the type of collision based on the cross sections
%% take into account energy loss in inelastic collisions
%% generate scattering and azimuth angles
%% in case of ionization handle the 'new' electron
t0 = sigma_e(1, energy_index);
t1 = t0 + sigma_e(2, energy_index);
t2 = t1 + sigma_e(3, energy_index);
rnd = rand;

if rnd < (t0 / t2) % elastic scattering
    chi = acos(1.0 - 2.0 * rand); % isotropic scattering
    eta = 2 * pi()*rand; % azimuthal angle
elseif rnd < (t1 / t2) % excitation
    energy = 0.5 * E_MASS * g * g; % electron energy
    energy = abs(energy - E_EXC_TH * EV_TO_J); % subtract energy loss for excitation
    g = sqrt(2.0 * energy / E_MASS); % relative velocity after energy loss
    chi = acos(1.0 - 2.0 * rand); % isotropic scattering
    eta = 2*pi() * rand; % azimuthal angle
else % ionization
    energy = 0.5 * E_MASS * g * g; % electron energy
    energy = abs(energy - E_ION_TH * EV_TO_J); % subtract energy loss for ionization
    e_ej = 10.0 * tan(rand * atan(energy / EV_TO_J / 20.0)) * EV_TO_J; % energy of the emitted electron
    e_sc = abs(energy - e_ej); % energy of incoming electron after collision
    g = sqrt(2.0 * e_sc / E_MASS); % relative velocity of incoming (original) electron
    g2 = sqrt(2.0 * e_ej / E_MASS); % relative velocity of emitted (new) electron
    chi = acos(sqrt(e_sc / energy)); % scattering angle for incoming electron
    chi2 = acos(sqrt(e_ej / energy)); % scattering angle for emitted electrons
    eta = 2*pi() * rand; % azimuthal angle for incoming electron
    eta2 = eta + pi; % azimuthal angle for emitted electron
    sc = sin(chi2);
    cc = cos(chi2);
    se = sin(eta2);
    ce = cos(eta2);
    gx = g2 * (ct * cc - st * sc * ce);
    gy = g2 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g2 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    N_e = N_e + 1; % add new electron
    WEIGHT_e(end+1) = WEIGHT_e(k+1);
    x_e(end+1) = xe;
    vx_e(end+1) = wx + F2 * gx;
    vy_e(end+1) = wy + F2 * gy;
    vz_e(end+1) = wz + F2 * gz;
    N_i = N_i + 1; % add new ion
    WEIGHT_i(end+1) = WEIGHT_e(k+1);
    x_i(end+1) = xe;
    vx_i(end+1) = RMB.random() ; % velocity is sampled from background thermal distribution
    vy_i(end+1) = RMB.random() ;
    vz_i(end+1) = RMB.random() ;
end
%% scatter the primary electron
sc = sin(chi);
cc = cos(chi);
se = sin(eta);
ce = cos(eta);
%% compute new relative velocity:
gx = g * (ct * cc - st * sc * ce);
gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

%% post-collision velocity of the colliding electron
vxe = wx + F2 * gx;
vye = wy + F2 * gy;
vze = wz + F2 * gz;
x_e(k+1)=xe; vx_e(k+1)=vxe; vy_e(k+1)=vye; vz_e(k+1)=vze; 

