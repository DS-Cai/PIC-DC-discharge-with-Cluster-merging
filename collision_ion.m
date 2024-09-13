vx_1=vx_i(k+1); vy_1=vy_i(k+1); vz_1=vz_i(k+1); vx_2=vx_a; vy_2=vy_a; vz_2=vz_a; 
%% calculate relative velocity before collision
%% random Maxwellian target atom already selected (vx_2,vy_2,vz_2 velocity components of target atom come with the call)
gx = vx_1 - vx_2;
gy = vy_1 - vy_2;
gz = vz_1 - vz_2;
g = sqrt(gx^2 + gy^2 + gz^2);
wx = 0.5 * (vx_1 + vx_2);
wy = 0.5 * (vy_1 + vy_2);
wz = 0.5 * (vz_1 + vz_2);
if gx == 0
    theta = 0.5 * pi;
else
    theta = atan2(sqrt(gy^2 + gz^2), gx);
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
t1 = sigma_i(1, energy_index);
t2 = t1 + sigma_i(2, energy_index);
rnd = rand;
if rnd < (t1 / t2)
    chi = acos(1.0 - 2.0 * rand);
else
    chi = pi;
end
eta = 2.0 * pi * rand;
sc = sin(chi);
cc = cos(chi);
se = sin(eta);
ce = cos(eta);
st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);
gx = g * (ct * cc - st * sc * ce);
gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
vx_1 = wx + 0.5 * gx;
vy_1 = wy + 0.5 * gy;
vz_1 = wz + 0.5 * gz;
vx_i(k+1)=vx_1; vy_i(k+1)=vy_1; vz_i(k+1)=vz_1;
