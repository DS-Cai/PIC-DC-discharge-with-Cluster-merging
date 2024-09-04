% 打开文件
f = fopen('OUTPUT/info.txt', 'w');
line_= repmat('-', 1, 80);
% 设置输出格式
formatSpec = '%.4e';
WEIGHT=1e5;
density = cumul_e_density(N_G / 2+1) / floor(no_of_cycles) /floor(N_T);  % e density @ center
plas_freq = E_CHARGE * sqrt(density / EPSILON0 / E_MASS);  % e plasma frequency @ center
meane = double(mean_energy_accu_center / floor(mean_energy_counter_center));  % e mean energy @ center
kT = 2.0 * meane * EV_TO_J / 3.0;  % k T_e @ center (approximate)
debye_length = sqrt(double(EPSILON0 * 2.0 * meane * EV_TO_J / 3.0 / density)) / E_CHARGE;  % e Debye length @ center
sim_time = no_of_cycles / FREQUENCY;  % simulated time
ecoll_freq = N_e_coll / sim_time / N_e;  % e collision frequency
icoll_freq = N_i_coll / sim_time / N_i;  % ion collision frequency

% 输出到文件
fprintf(f, '########################## eduPIC simulation report ############################\n');
fprintf(f, 'Simulation parameters:\n');
fprintf(f, 'Gap distance                          = %f [m]\n', L);
fprintf(f, '# of grid divisions                   = %d\n', N_G);
fprintf(f, 'Frequency                             = %f [Hz]\n', FREQUENCY);
fprintf(f, '# of time steps / period              = %d\n', N_T);
fprintf(f, '# of electron / ion time steps        = %d\n', N_SUB);
fprintf(f, 'Voltage amplitude                     = %f [V]\n', VOLTAGE);
fprintf(f, 'Pressure (Ar)                         = %f [Pa]\n', PRESSURE);
fprintf(f, 'Temperature                           = %f [K]\n', TEMPERATURE);
fprintf(f, 'Superparticle weight                  = %f\n', WEIGHT);
fprintf(f, '# of simulation cycles in this run    = %d\n', no_of_cycles);
fprintf(f, '%s\n', line_);
fprintf(f, 'Plasma characteristics:\n');
fprintf(f, 'Electron density @ center             = %e [m^{-3}]\n', density);
fprintf(f, 'Plasma frequency @ center             = %e [rad/s]\n', plas_freq);
fprintf(f, 'Debye length @ center                 = %e [m]\n', debye_length);
fprintf(f, 'Electron collision frequency          = %e [1/s]\n', ecoll_freq);
fprintf(f, 'Ion collision frequency               = %e [1/s]\n', icoll_freq);
fprintf(f, '%s\n', line_);
fprintf(f, 'Stability and accuracy conditions:\n');
conditions_OK = true;
c = plas_freq * DT_E;
fprintf(f, 'Plasma frequency @ center * DT_e      = %e (OK if less than 0.20)\n', c);
if (c > 0.2)
    conditions_OK = false;
end
c = DX / debye_length;
fprintf(f, 'DX / Debye length @ center            = %e (OK if less than 1.00)\n', c);
if (c > 1.0)
    conditions_OK = false;
end
c = max_electron_coll_freq(CS_RANGES,DE_CS,EV_TO_J,E_MASS,sigma_tot_e)* DT_E;
fprintf(f, 'Max. electron coll. frequency * DT_E  = %e (OK if less than 0.05)\n', c);
if (c > 0.05)
    conditions_OK = false;
end
c = max_ion_coll_freq(CS_RANGES,DE_CS,EV_TO_J,MU_ARAR,sigma_tot_i) * DT_I;
fprintf(f, 'Max. ion coll. frequency * DT_I       = %e (OK if less than 0.05)\n', c);
if (c > 0.05)
    conditions_OK = false;
end
if conditions_OK == false
    fprintf(f, line_);
    fprintf(f, '** STABILITY AND ACCURACY CONDITION(S) VIOLATED - REFINE SIMULATION SETTINGS! **\n');
    fprintf(f, '%s\n', line_);
    fclose(f);
    error('** STABILITY AND ACCURACY CONDITION(S) VIOLATED - REFINE SIMULATION SETTINGS! **\n')
    fprintf(f, '>> eduPIC: ERROR: STABILITY AND ACCURACY CONDITION(S) VIOLATED! \n');
    fprintf(f, '>> eduPIC: for details see ''info.txt'' and refine simulation settings!\n');
else
    % calculate maximum energy for which the Courant condition holds:
    v_max = DX / DT_E;
    e_max = 0.5 * E_MASS * v_max^2 / EV_TO_J;
    fprintf(f, 'Max e- energy for CFL     condition   = %f\n', e_max);
    fprintf(f, 'Check EEPF to ensure that CFL is fulfilled for the majority of the electrons!\n');
    fprintf(f, line_);

    % saving of the following data is done here as some of the further lines need data
    % that are computed / normalized in these functions

    fprintf('>> eduPIC: saving diagnostics data\n');
    %% save_density;
    ff = fopen('OUTPUT/density.dat', 'w');
    fprintf(ff, '%.12e\n', cumul_e_density * c);
    fprintf(ff, '%.12e\n', cumul_i_density * c);
    fclose(ff);
    %% save_eepf
    ff = fopen('OUTPUT/eepf.dat', 'w');
    h = sum(eepf) * DE_EEPF;
    fprintf(ff, '%.12e\n', [((0:N_EEPF-1) + 0.5) * DE_EEPF; eepf./(h.*sqrt(((0:N_EEPF-1) + 0.5) * DE_EEPF))]);
    fclose(ff);
    %% save_ifed
    ff = fopen('OUTPUT/ifed.dat', 'w');
    h_pow = sum(ifed_pow) * DE_IFED;
    h_gnd = sum(ifed_gnd) * DE_IFED;
    fprintf(ff, '%.12e\n', [((0:N_IFED-1) + 0.5) * DE_IFED; ifed_pow/h_pow; ifed_gnd/h_gnd]);
    fclose(ff);
    mean_i_energy_pow = sum(((0:N_IFED-1) + 0.5) * DE_IFED .* (ifed_pow/h_pow));
    mean_i_energy_gnd = sum(((0:N_IFED-1) + 0.5) * DE_IFED .* (ifed_gnd/h_gnd));
    %% norm_all_xt();
    f1 = double(N_XT) / double(no_of_cycles * N_T);
    f2 = WEIGHT / (ELECTRODE_AREA * DX) / (no_of_cycles * (PERIOD / double(N_XT)));
    pot_xt = f1 * pot_xt;
    efield_xt = f1 * efield_xt;
    ne_xt = f1 * ne_xt;
    ni_xt = f1 * ni_xt;
    ue_xt = arrayfun(@(x, y) (y > 0) * (x / y), ue_xt, counter_e_xt);
    je_xt = -ue_xt .* ne_xt * E_CHARGE;
    meanee_xt = arrayfun(@(x, y) (y > 0) * (x / y), meanee_xt, counter_e_xt);
    ioniz_rate_xt = arrayfun(@(x, y) (y > 0) * (x * f2), ioniz_rate_xt, counter_e_xt);
    ui_xt = arrayfun(@(x, y) (y > 0) * (x / y), ui_xt, counter_i_xt);
    ji_xt = ui_xt .* ni_xt * E_CHARGE;
    meanei_xt = arrayfun(@(x, y) (y > 0) * (x / y), meanei_xt, counter_i_xt);
    powere_xt = je_xt .* efield_xt;
    poweri_xt = ji_xt .* efield_xt;
    %% save_all_xt
    save_xt_1(pot_xt, 'OUTPUT\pot_xt.dat');
    save_xt_1(efield_xt, 'OUTPUT\efield_xt.dat');
    save_xt_1(ne_xt, 'OUTPUT\ne_xt.dat');
    save_xt_1(ni_xt, 'OUTPUT\ni_xt.dat');
    save_xt_1(je_xt, 'OUTPUT\je_xt.dat');
    save_xt_1(ji_xt, 'OUTPUT\ji_xt.dat');
    save_xt_1(powere_xt, 'OUTPUT\powere_xt.dat');
    save_xt_1(poweri_xt, 'OUTPUT\poweri_xt.dat');
    save_xt_1(meanee_xt, 'OUTPUT\meanee_xt.dat');
    save_xt_1(meanei_xt, 'OUTPUT\meanei_xt.dat');
    save_xt_1(ioniz_rate_xt, 'OUTPUT\ioniz_xt.dat');
    %% 
    fprintf(f, 'Particle characteristics at the electrodes:\n');
    fprintf(f, 'Ion flux at powered electrode         = %e [m^{-2} s^{-1}]\n', N_i_abs_pow * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD));
    fprintf(f, 'Ion flux at grounded electrode        = %e [m^{-2} s^{-1}]\n', N_i_abs_gnd * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD));
    fprintf(f, 'Mean ion energy at powered electrode  = %f [eV]\n', mean_i_energy_pow);
    fprintf(f, 'Mean ion energy at grounded electrode = %f [eV]\n', mean_i_energy_gnd);
    fprintf(f, 'Electron flux at powered electrode    = %e [m^{-2} s^{-1}]\n', N_e_abs_pow * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD));
    fprintf(f, 'Electron flux at grounded electrode   = %e [m^{-2} s^{-1}]\n', N_e_abs_gnd * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD));

    % calculate spatially and temporally averaged power absorption by the electrons and ions

    power_e = sum(powere_xt(:)) / (N_XT * N_G);
    power_i = sum(poweri_xt(:)) / (N_XT * N_G);
    fprintf(f, line_);
    fprintf(f, 'Absorbed power calculated as <j*E>:\n');
    fprintf(f, 'Electron power density (average)      = %f [W m^{-3}]\n', power_e);
    fprintf(f, 'Ion power density (average)           = %f [W m^{-3}]\n', power_i);
    fprintf(f, 'Total power density (average)         = %f [W m^{-3}]\n', power_e + power_i);
    fprintf(f, line_);
    fclose(f);
end