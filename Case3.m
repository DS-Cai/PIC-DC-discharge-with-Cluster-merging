clear all; clc; fclose all;
parameter_define;
measurement_mode=1; %Record the variation of variables with time
sigma_e=set_electron_cross_sections_ar(CS_RANGES,DE_CS,E_EXC_TH,E_ION_TH);
sigma_i=set_ion_cross_sections_ar(CS_RANGES,DE_CS);
sigma_tot_e=(sigma_e(1, :) + sigma_e(2, :) + sigma_e(3, :)) .* GAS_DENSITY;
sigma_tot_i = (sigma_i(1, :) + sigma_i(2, :)) .* GAS_DENSITY; 
%% Initialization
no_of_cycles=1;
cycle = 1; 
merge=1; % If merge
test_minsize=[8 12 16 20];
min_size=test_minsize(1); %
Init; %Initialize
fprintf("Running initializing cycle\n");
Time = 0;
%% Main cycle
do_one_cycle();
save_particle_data;
check_and_save_info;
fprintf('Simulation of %d cycle(s) is completed.\n', no_of_cycles);
OUTPUT;
