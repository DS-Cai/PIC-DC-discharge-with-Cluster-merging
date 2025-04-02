function save_xt_1(distr, fname)
    f = fopen(fname, 'w');
    fprintf(f, '%.8e ', distr(:));
    fprintf(f, '\n');
    fclose(f);
end