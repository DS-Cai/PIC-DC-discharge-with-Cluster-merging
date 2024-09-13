%% constants
PI = 3.141592653589793;      % mathematical constant Pi
TWO_PI = 2.0 * PI;               % two times Pi
E_CHARGE = 1.60217662e-19;         % electron charge [C]
EV_TO_J = E_CHARGE;               % eV <-> Joule conversion factor
E_MASS = 9.10938356e-31;         % mass of electron [kg]
AR_MASS = 6.63352090e-26;         % mass of argon atom [kg]
MU_ARAR = AR_MASS / 2.0;          % reduced mass of two argon atoms [kg]
K_BOLTZMANN = 1.38064852e-23;         % Boltzmann's constant [J/K]
EPSILON0 = 8.85418781e-12;         % permittivity of free space [F/m]
%% simulation parameters
N_G = 1000;                    % number of grid points
N_T = 2e4;                   % time steps within an RF period
FREQUENCY = 13.56e6;                % driving frequency [Hz]
VOLTAGE = 1000;                  % voltage amplitude [V]
L = 0.01;                  % electrode gap [m]
PRESSURE = 133.33;                   % gas pressure [Pa]
TEMPERATURE = 350.0;                  % background gas temperature [K]
ELECTRODE_AREA = 1e-4;                 % (fictive) electrode area [m^2]
N_INIT = 10000;                   % number of initial electrons and ions
%% additional (derived) constants 
PERIOD = 1.0 / FREQUENCY;                           % RF period length [s]
DT_E = 5e-13;         % electron time step [s]
N_SUB = 20;                                        % ions move only in these cycles (subcycling)
DT_I = N_SUB * DT_E;                              % ion time step [s]
DX = L /(N_G - 1);          % spatial grid division [m]
INV_DX = 1.0 / DX;                                  % inverse of spatial grid size [1/m]
GAS_DENSITY = PRESSURE / (K_BOLTZMANN * TEMPERATURE);    % background gas density [1/m^3]
OMEGA = TWO_PI * FREQUENCY;                        % angular frequency [rad/s]
%% electron and ion cross sections
N_CS = 5;                      % total number of processes / cross sections
E_ELA = 0;                      % process identifier: electron/elastic
E_EXC = 1;                      % process identifier: electron/excitation
E_ION = 2;                      % process identifier: electron/ionization
I_ISO = 3;                      % process identifier: ion/elastic/isotropic
I_BACK = 4;                      % process identifier: ion/elastic/backscattering
E_EXC_TH = 11.5;                   % electron impact excitation threshold [eV]
E_ION_TH = 15.8;                   % electron impact ionization threshold [eV]
CS_RANGES = 1000000;                % number of entries in cross section arrays
DE_CS = 0.001;                  % energy division in cross section arrays [eV]
sigma=zeros(1,CS_RANGES);       %set of cross section arrays
sigma_tot_e=zeros(1,CS_RANGES);       %total macroscopic cross section of electrons
sigma_tot_i=zeros(1,CS_RANGES);       %total macroscopic cross section of ions
%% particle coordinates
N_e = 0; % number of electrons
N_i = 0; % number of ions
x_e=0; vx_e=0; vy_e=0; vz_e=0; %coordinates of electrons (one spatial, three velocity components)
x_i=0; vx_i=0; vy_i=0; vz_i=0; %coordinates of ions (one spatial, three velocity components) 
efield=zeros(1, N_G);  % electric field and potential
pot=zeros(1, N_G);     % electric field and potential
e_density=zeros(1, N_G);  % electron and ion densities
i_density=zeros(1, N_G); 
xN_e=zeros(1, N_G); 
xN_i=zeros(1, N_G); 
xG_e=zeros(N_G,1); 
xG_i=zeros(N_G,1); 
cumul_e_density=zeros(1, N_G); %cumulative densities
cumul_i_density=zeros(1, N_G);
N_e_abs_pow = uint64(0); % counter for electrons absorbed at the powered electrode
N_e_abs_gnd = uint64(0); % counter for electrons absorbed at the grounded electrode
N_i_abs_pow = uint64(0); % counter for ions absorbed at the powered electrode
N_i_abs_gnd = uint64(0); % counter for ions absorbed at the grounded electrode
%% electron energy probability function
N_EEPF = 2000;   % number of energy bins in Electron Energy Probability Function (EEPF)
DE_EEPF = 0.05;  % resolution of EEPF [eV]
eepf=zeros(1,N_EEPF);  % time integrated EEPF in the center of the plasma
%% ion flux-energy distributions
N_IFED = 100;  % number of energy bins in Ion Flux-Energy Distributions (IFEDs)
DE_IFED = 1.0; % resolution of IFEDs [eV]
ifed_pow=zeros(1,N_IFED);  % IFED at the powered electrode
ifed_gnd=zeros(1,N_IFED);  % IFED at the grounded electrode
mean_i_energy_pow=0;  % mean ion energy at the powered electrode
mean_i_energy_gnd=0;  % mean ion energy at the grounded electrode
%% spatio-temporal (XT) distributions
N_BIN = 20;   % number of time steps binned for the XT distributions
N_XT = N_T / N_BIN;  % number of spatial bins for the XT distributions
pot_xt = zeros(1, N_G* N_XT);             % XT distribution of the potential
efield_xt = zeros(1, N_G* N_XT);             % XT distribution of the electric field
ne_xt = zeros(1, N_G* N_XT);             % XT distribution of the electron density
ni_xt = zeros(1, N_G* N_XT);             % XT distribution of the ion density
ue_xt = zeros(1, N_G* N_XT);             % XT distribution of the mean electron velocity
ui_xt = zeros(1, N_G* N_XT);             % XT distribution of the mean ion velocity
je_xt = zeros(1, N_G* N_XT);             % XT distribution of the electron current density
ji_xt = zeros(1, N_G* N_XT);             % XT distribution of the ion current density
powere_xt = zeros(1, N_G* N_XT);             % XT distribution of the electron powering (power absorption) rate
poweri_xt = zeros(1, N_G* N_XT);             % XT distribution of the ion powering (power absorption) rate
meanee_xt = zeros(1, N_G* N_XT);             % XT distribution of the mean electron energy
meanei_xt = zeros(1, N_G* N_XT);             % XT distribution of the mean ion energy
counter_e_xt = zeros(1, N_G* N_XT);             % XT counter for electron properties
counter_i_xt = zeros(1, N_G* N_XT);             % XT counter for ion properties
ioniz_rate_xt = zeros(1, N_G* N_XT);             % XT distribution of the ionisation rate
mean_energy_accu_center=0;                     % mean electron energy accumulator in the center of the gap               
mean_energy_counter_center = uint64(0);                 % mean electron energy counter in the center of the gap
N_e_coll = uint64(0);                 % counter for electron collisions
N_i_coll = uint64(0);                 % counter for ion collisions
Time=0;  % total simulated time (from the beginning of the simulation)
cycles_done=0; no_of_cycles=0; % current cycle and total cycles in the run, cycles completed (from the beginning of the simulation)
arg1=0;  % used for reading command line arguments
measurement_mode=0; % flag that controls measurements and data saving
 %%
 % 创建正态分布对象
 rng('default');  % 设置随机数种子为默认值
 MTgen = rng('shuffle', 'twister');
 RMB = makedist('Normal', 0, sqrt(K_BOLTZMANN * TEMPERATURE / AR_MASS));
