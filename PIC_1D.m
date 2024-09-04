clear all; clc; fclose all;
parameter_define;
measurement_mode = 1; % Save the time-dependent variables
sigma_e=set_electron_cross_sections_ar(CS_RANGES,DE_CS,E_EXC_TH,E_ION_TH);
sigma_i=set_ion_cross_sections_ar(CS_RANGES,DE_CS);
sigma_tot_e=(sigma_e(1, :) + sigma_e(2, :) + sigma_e(3, :)) .* GAS_DENSITY;
sigma_tot_i = (sigma_i(1, :) + sigma_i(2, :)) .* GAS_DENSITY; 
%% Initialize
no_of_cycles=1;
cycle = 1; 
Init; %Initialize
Time = 0;
%% Cycle
do_one_cycle();
save_particle_data;
if measurement_mode
    check_and_save_info;
end
fprintf('eduPIC: simulation of %d cycle(s) is completed.\n', no_of_cycles);
