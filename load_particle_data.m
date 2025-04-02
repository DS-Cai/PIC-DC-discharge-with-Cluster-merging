% load particle coordinates  
fid = fopen('OUTPUT/picdata.bin', 'rb');
if fid == -1
    error('ERROR: No particle data file found, try running initial cycle using argument')
else
Time = fread(fid, 1, 'double');
N_e = fread(fid, 1, 'int');
WEIGHT_e = fread(fid, N_e, 'double');
x_e = fread(fid, N_e, 'double');
vx_e = fread(fid, N_e, 'double');
vy_e = fread(fid, N_e, 'double');
vz_e = fread(fid, N_e, 'double');
N_i = fread(fid, 1, 'int');
WEIGHT_i = fread(fid, N_i, 'double');
x_i = fread(fid, N_i, 'double');
vx_i = fread(fid, N_i, 'double');
vy_i = fread(fid, N_i, 'double');
vz_i = fread(fid, N_i, 'double');

fclose(fid);

disp(['>> eduPIC: data loaded : ', num2str(N_e), ' electrons ', num2str(N_i), ' ions, ', ' time before is ', num2str(Time), ' [s]']);
end
