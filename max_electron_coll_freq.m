function nu_max = max_electron_coll_freq(CS_RANGES,DE_CS,EV_TO_J,E_MASS,sigma_tot_e)
    nu_max = 0;
    for i = 1:CS_RANGES
        e = (i - 1) * DE_CS;
        v = sqrt(2.0 * e * EV_TO_J / E_MASS);
        nu = v * sigma_tot_e(i);
        if nu > nu_max
            nu_max = nu;
        end
    end
end