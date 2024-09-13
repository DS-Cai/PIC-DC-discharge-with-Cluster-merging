f = fopen('OUTPUT\picdata.bin', 'wb');
fwrite(f, Time, 'double');
fwrite(f, N_e, 'int');
fwrite(f, WEIGHT_e, 'double');
fwrite(f, x_e, 'double');
fwrite(f, vx_e, 'double');
fwrite(f, vy_e, 'double');
fwrite(f, vz_e, 'double');
fwrite(f, N_i, 'int');
fwrite(f, WEIGHT_i, 'double');
fwrite(f, x_i, 'double');
fwrite(f, vx_i, 'double');
fwrite(f, vy_i, 'double');
fwrite(f, vz_i, 'double');
fclose(f);

disp(['>>Data saved : ', num2str(N_e), ' electrons ', num2str(N_i), ' ions, ', ' time is ', num2str(Time, '%.2e'), ' [s]']);
