function nu_max = max_ion_coll_freq(CS_RANGES,DE_CS,EV_TO_J,MU_ARAR,sigma_tot_i)
    nu_max = 0;
    for i = 1:CS_RANGES
        e = (i - 1) * DE_CS;
        g = sqrt(2.0 * e * EV_TO_J / MU_ARAR);
        nu = g * sigma_tot_i(i);
        if nu > nu_max
            nu_max = nu;
        end
    end
end