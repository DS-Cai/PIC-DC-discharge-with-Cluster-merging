function sigma=set_electron_cross_sections_ar(CS_RANGES,DE_CS,E_EXC_TH,E_ION_TH)
    disp("eduPIC: Setting e- / Ar cross sections");
    e = zeros(1, CS_RANGES);
    sigma = zeros(3, numel(e));
for i=1: numel(e)
    if i==1
        e(i)=DE_CS;
    else
    e(i) = DE_CS * (i-1);
    end
    sigma(1, i) =1e-20 * (abs(6.0 ./ (1.0 + (e(i) / 0.1) + (e(i) / 0.6).^2).^3.3 ...
        - 1.1 * e(i).^1.4 ./ (1.0 + (e(i) / 15.0).^1.2) ./ sqrt(1.0 + (e(i) / 5.5).^2.5 + (e(i) / 60.0).^4.1)) ...
        + 0.05 ./ (1.0 + e(i) / 10.0).^2.0 + 0.01 * e(i).^3.0 ./ (1.0 + (e(i) / 12.0).^6.0));
    if e(i)> E_EXC_TH
    sigma(2, i) =1e-20 * (0.034 * (e(i) - 11.5).^1.1 .* (1.0 + (e(i) / 15.0).^2.8) ./ (1.0 + (e(i) / 23.0).^5.5) ...
            + 0.023 * (e(i) - 11.5) ./ (1.0 + e(i) / 80.0).^1.9);
    else
        sigma(2, i)=0;
    end
    if e(i)> E_ION_TH
    sigma(3, i) = 1e-20 * (970.0 * (e(i) - 15.8) ./ (70.0 + e(i)).^2.0+ 0.06 * (e(i) - 15.8).^2.0 .* exp(-e(i) / 9));
    else
        sigma(3, i)=0;
    end
end