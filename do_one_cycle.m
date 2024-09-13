%% Initialize variables
DV = ELECTRODE_AREA * DX;
FACTOR_E = DT_E / E_MASS * E_CHARGE;
FACTOR_I = DT_I / AR_MASS * E_CHARGE;
MIN_X = 0.45 * L; MAX_X = 0.55 * L; % min. /max. position for EEPF collection
k = 0; t = 0; p = 0; energy_index = 0; rmod = 0; rint = 0; g = 0; g_sqr = 0;
gx = 0; gy = 0; gz = 0; vx_a = 0; vy_a = 0; vz_a = 0; e_x = 0; energy = 0;
nu = 0; p_coll = 0; v_sqr = 0; velocity = 0; mean_v = 0; rate = 0; out = false;
rho = []; t_index = 0; CPU_t=zeros(N_T+1,1);
%% 
for t = 0:N_T-1
    tic
    Time = Time + DT_E;  % update of the total simulated time
    t_index = floor(t/ N_BIN);  % index for XT distributions
%% step 1: compute densities at grid points
    e_density(:) =0;  % electron density - computed in every time step
    for k = 0:N_e-1  % calculate electron density on cell nodes
        FACTOR_W = WEIGHT_e(k+1) / DV;
        rmod = mod(x_e(k+1) * INV_DX, 1); % Fractional Part. get superelectron coordinates
        p = floor(x_e(k+1) * INV_DX); % Integer Part
        e_density(p+1) = e_density(p+1) + (1 - rmod) * FACTOR_W;
        e_density(p+2) = e_density(p+2) + rmod * FACTOR_W;
    end
    e_density(1) = e_density(1) * 2;
    e_density(end) = e_density(end) * 2;
    cumul_e_density = cumul_e_density + e_density;
    if mod(t, N_SUB) == 0
        i_density(:) = 0;     % ion density - computed in every N_SUB-th time steps (subcycling)
        for k = 0:N_i-1
            FACTOR_W = WEIGHT_i(k+1) / DV;
            rmod = mod(x_i(k+1) * INV_DX, 1);
            p = fix(x_i(k+1) * INV_DX);
            i_density(p+1) = i_density(p+1) + (1 - rmod) * FACTOR_W;
            i_density(p+2) = i_density(p+2) + rmod * FACTOR_W;
        end
        i_density(1) = i_density(1) * 2.0;
        i_density(end) = i_density(end) * 2.0;
    end
    cumul_i_density = cumul_i_density + i_density;   
%% step 2: solve Poisson equation   
% get charge density
    rho = E_CHARGE * (i_density - e_density);
    [efield,pot]=solve_Poisson(rho,Time,DX,EPSILON0,N_G,VOLTAGE,OMEGA,INV_DX);
%% steps 3 & 4: move particles according to electric field interpolated to particle positions
for k = 0:N_e-1                      % move all electrons in every time step
    rmod = mod(x_e(k+1) * INV_DX, 1);
    p = floor(x_e(k+1) * INV_DX);
    e_x = (1.0 - rmod) * efield(p+1) + rmod* efield(p+2);
    if measurement_mode
        % measurements: 'x' and 'v' are needed at the same time, i.e. old
        % 'x' and mean 'v' 
        mean_v = vx_e(k+1) - 0.5 * e_x * FACTOR_E;
        counter_e_xt(p * N_XT + t_index+1) = counter_e_xt(p * N_XT + t_index+1) +  (1.0 - rmod)*WEIGHT_e(k+1);
        counter_e_xt((p + 1) * N_XT + t_index+1) = counter_e_xt((p + 1) * N_XT + t_index+1) + rmod*WEIGHT_e(k+1);
        ue_xt(p * N_XT + t_index+1) = ue_xt(p * N_XT + t_index+1) + (1.0 - rmod) * mean_v;
        ue_xt((p + 1) * N_XT + t_index+1) = ue_xt((p + 1) * N_XT + t_index+1) + rmod * mean_v;
        v_sqr = mean_v * mean_v + vy_e(k+1) * vy_e(k+1) + vz_e(k+1) * vz_e(k+1);
        energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
        meanee_xt(p * N_XT + t_index+1) = meanee_xt(p * N_XT + t_index+1) + (1.0 - rmod) * energy* WEIGHT_e(k+1);
        meanee_xt((p + 1) * N_XT + t_index+1) = meanee_xt((p + 1) * N_XT + t_index+1) + rmod * energy* WEIGHT_e(k+1);
        energy_index = min( floor(energy / DE_CS), CS_RANGES - 1)+1;
        velocity = sqrt(v_sqr);
        rate = sigma_e(3, energy_index) * velocity * DT_E * GAS_DENSITY;
        ioniz_rate_xt(p * N_XT + t_index+1) = ioniz_rate_xt(p * N_XT + t_index+1) + (1.0 - rmod) * rate;
        ioniz_rate_xt((p + 1) * N_XT + t_index+1) = ioniz_rate_xt((p + 1) * N_XT + t_index+1) + rmod * rate;
        % measure EEPF in the center
            energy_index =  floor(energy / DE_EEPF)+1;
            if energy_index < N_EEPF+1
                eepf(energy_index) = eepf(energy_index) + 1.0*WEIGHT_e(k+1);
            end
            mean_energy_accu_center = mean_energy_accu_center + energy;
            mean_energy_counter_center = mean_energy_counter_center + 1;
    end
    % update velocity and position
    vx_e(k+1) = vx_e(k+1) - e_x * FACTOR_E;
    x_e(k+1) = x_e(k+1) + vx_e(k+1) * DT_E;
end

if mod(t, N_SUB) == 0     % move all ions in every N_SUB-th time steps (subcycling)
    for k = 0:N_i-1
        rmod = mod(x_i(k+1) * INV_DX, 1);
        p = fix(x_i(k+1) * INV_DX);
        e_x = (1.0 - rmod) * efield(p+1) + rmod * efield(p + 2);
        if measurement_mode
            % measurements: 'x' and 'v' are needed at the same time, i.e. old 'x' and mean 'v'
            mean_v = vx_i(k+1) + 0.5 * e_x * FACTOR_I;
            counter_i_xt(p * N_XT + t_index+1) = counter_i_xt(p * N_XT + t_index+1) + (1.0 - rmod)* WEIGHT_i(k+1);
            counter_i_xt((p + 1) * N_XT + t_index+1) = counter_i_xt((p + 1) * N_XT + t_index+1) + rmod* WEIGHT_i(k+1);
            ui_xt(p * N_XT + t_index+1) = ui_xt(p * N_XT + t_index+1) + (1.0 - rmod) * mean_v;
            ui_xt((p + 1) * N_XT + t_index+1) = ui_xt((p + 1) * N_XT + t_index+1) + rmod * mean_v;
            v_sqr = mean_v * mean_v + vy_i(k+1) * vy_i(k+1) + vz_i(k+1) * vz_i(k+1);
            energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
            meanei_xt(p * N_XT + t_index+1) = meanei_xt(p * N_XT + t_index+1) + (1.0 - rmod) * energy* WEIGHT_i(k+1);
            meanei_xt((p + 1) * N_XT + t_index+1) = meanei_xt((p + 1) * N_XT + t_index+1) + rmod * energy* WEIGHT_i(k+1);
        end
        % update velocity and position
        vx_i(k+1) = vx_i(k+1) + e_x * FACTOR_I;
        x_i(k+1) = x_i(k+1) + vx_i(k+1) * DT_I;
    end
end
%% step 5: check boundaries
k = 0;
while k < length(x_e)
    out = false;
    if x_e(k+1) < 0
        N_e_abs_pow = N_e_abs_pow + 1;
        out = true; % the electron is out at the powered electrode
    end
    if x_e(k+1) > L
        N_e_abs_gnd = N_e_abs_gnd + 1;
        out = true; % the electron is out at the grounded electrode
    end
    if out
        % remove the electron if out
        WEIGHT_e(k+1)=WEIGHT_e(end);
        WEIGHT_e(end)=[];
        x_e(k+1) = x_e(end);
        x_e(end) = [];
        vx_e(k+1) = vx_e(end);
        vx_e(end) = [];
        vy_e(k+1) = vy_e(end);
        vy_e(end) = [];
        vz_e(k+1) = vz_e(end);
        vz_e(end) = [];
        N_e = N_e - 1;
    else
        k = k + 1;
    end
end

if mod(t, N_SUB) == 0      % check boundaries for all ions in every N_SUB-th time steps (subcycling)
k = 0;
while k < length(x_i)
    out = false;
    if x_i(k+1) < 0
        N_i_abs_pow = N_i_abs_pow + 1;
        out = true; % the ion is out at the powered electrode
        v_sqr = vx_i(k+1) * vx_i(k+1) + vy_i(k+1) * vy_i(k+1) + vz_i(k+1) * vz_i(k+1);
        energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
        energy_index = floor(energy / DE_IFED)+1;
        if energy_index < N_IFED+1
            ifed_pow(energy_index) = ifed_pow(energy_index) + 1*WEIGHT_i(k+1); % save IFED at the powered electrode
        end
    end
    if x_i(k+1) > L
        N_i_abs_gnd = N_i_abs_gnd + 1;
        out = true; % the ion is out at the grounded electrode
        v_sqr = vx_i(k+1) * vx_i(k+1) + vy_i(k+1) * vy_i(k+1) + vz_i(k+1) * vz_i(k+1);
        energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
        energy_index = floor(energy / DE_IFED)+1;
        if energy_index < N_IFED+1
            ifed_gnd(energy_index) = ifed_gnd(energy_index) + 1; % save IFED at the grounded electrode
        end
    end
    if out
        % delete the ion if out
        WEIGHT_i(k+1)=WEIGHT_i(end);
        WEIGHT_i(end)=[];
        x_i(k+1) = x_i(end);
        x_i(end) = [];
        vx_i(k+1) = vx_i(end);
        vx_i(end) = [];
        vy_i(k+1) = vy_i(end);
        vy_i(end) = [];
        vz_i(k+1) = vz_i(end);
        vz_i(end) = [];
        N_i = N_i - 1;
    else
        k = k + 1;
    end
end
end
%% step 6: collisions
for k = 0:N_e-1
    p = floor(x_e(k+1) * INV_DX);
    v_sqr = vx_e(k+1) * vx_e(k+1) + vy_e(k+1) * vy_e(k+1) + vz_e(k+1) * vz_e(k+1);
    velocity = sqrt(v_sqr);
    energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
    energy_index = min(floor(energy / DE_CS + 0.5), CS_RANGES - 1)+1;
    nu = sigma_tot_e(energy_index) * velocity;
    p_coll = 1 - exp(-nu * DT_E); % collision probability for electrons
    if rand() < p_coll % electron collision takes place
        collision_electron
        N_e_coll = N_e_coll + 1;
    end
end
if mod(t, N_SUB) == 0
    for k = 0:N_i-1
        vx_a = random(RMB); % pick velocity components of a random gas atoms
        vy_a = random(RMB);
        vz_a = random(RMB);
        gx = vx_i(k+1) - vx_a; % compute the relative velocity of the collision partners
        gy = vy_i(k+1) - vy_a;
        gz = vz_i(k+1) - vz_a;
        g_sqr = gx * gx + gy * gy + gz * gz;
        g = sqrt(g_sqr);
        energy = 0.5 * MU_ARAR * g_sqr / EV_TO_J;
        energy_index = min(floor(energy / DE_CS + 0.5), CS_RANGES - 1)+1;
        nu = sigma_tot_i(energy_index) * g;
        p_coll = 1 - exp(-nu * DT_I); % collision probability for ions
        if rand() < p_coll % ion collision takes place
            collision_ion;
            N_i_coll = N_i_coll + 1;
        end
    end
end
%%
for p = 0:N_G-1
    pot_xt(p*N_XT + t_index + 1) = pot_xt(p*N_XT + t_index + 1) + pot(p+1);
    efield_xt(p*N_XT + t_index + 1) = efield_xt(p*N_XT + t_index + 1) + efield(p+1);
    ne_xt(p*N_XT + t_index + 1) = ne_xt(p*N_XT + t_index + 1) + e_density(p+1);
    ni_xt(p*N_XT + t_index + 1) = ni_xt(p*N_XT + t_index + 1) + i_density(p+1);
end
N_e_real=sum(WEIGHT_e)./1e5;
N_i_real=sum(WEIGHT_i)./1e5;
xG_e=zeros(N_G,1); 
xN_e=zeros(1, N_G); 
for k = 0:N_e-1
    p = floor(x_e(k+1) * INV_DX);
    xN_e(p+1)=xN_e(p+1)+1;
    xG_e(p+1,xN_e(p+1))=k+1;
end
if merge
[x_e,vx_e,vy_e,vz_e,WEIGHT_e]=merging(xN_e,xG_e,x_e,vx_e,vy_e,vz_e,WEIGHT_e,DX,N_G,8);
end
N_e=length(WEIGHT_e);
xG_i=zeros(N_G,1);
xN_i=zeros(1, N_G);
if mod(t, N_SUB) == 0
    for k = 0:N_i-1
        p = fix(x_i(k+1) * INV_DX);
        xN_i(p+1)=xN_i(p+1)+1;
        xG_i(p+1,xN_i(p+1))=k+1;
    end
    if merge
    [x_i,vx_i,vy_i,vz_i,WEIGHT_i]=merging(xN_i,xG_i,x_i,vx_i,vy_i,vz_i,WEIGHT_i,DX,N_G,8);
    end
    N_i=length(WEIGHT_i);
end
if mod(t, 100) == 0
    fprintf(' time = %8d t = %8d #e = %8d #e_real = %8d #i = %8d #i_real = %8d\n', Time, t, N_e,N_e_real, N_i,N_i_real);
end
datafile = fopen('OUTPUT/conv.dat', 'a');
fprintf(datafile, '%d %d %d\n', Time, N_e, N_i);
fclose(datafile);
CPU_t(t+2)=CPU_t(t+1)+toc;
end

    


