%% 初始化
N_e = N_INIT;
N_i =  N_INIT;
WEIGHT_e = 1e5.*ones(1, N_e);                  % weight of superparticles
WEIGHT_i = 1e5.*ones(1, N_i);                  % weight of superparticles
x_e = L * rand(1,N_INIT);
x_i = L * rand(1,N_INIT);
vx_e = zeros(1,N_INIT);
vy_e = zeros(1,N_INIT);
vz_e = zeros(1,N_INIT);
vx_i = zeros(1,N_INIT);
vy_i = zeros(1,N_INIT);
vz_i = zeros(1,N_INIT);
