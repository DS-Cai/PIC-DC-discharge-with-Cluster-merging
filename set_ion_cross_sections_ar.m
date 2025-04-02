function sigma=set_ion_cross_sections_ar(CS_RANGES,DE_CS)
    disp("Setting Ar+ / Ar cross sections");

    qiso = @(e_lab) 2e-19 * e_lab.^(-0.5) ./ (1.0 + e_lab) + ...
        3e-19 * e_lab ./ (1.0 + e_lab / 3.0).^2.0;

    qmom = @(e_lab) 1.15e-18 * e_lab.^(-0.1) .* (1.0 + 0.015 ./ e_lab).^0.6;

    qback = @(x) (qmom(x) - qiso(x)) / 2.0;

    e = zeros(1, CS_RANGES);
    e(1) = 2.0 * DE_CS;
    for i = 2:numel(e)
        e(i) = 2.0 * DE_CS * (i-1);
    end

    sigma = zeros(2, numel(e));

    for i = 1:numel(e)
        sigma(1, i) = qiso(e(i)); 
        sigma(2, i) = qback(e(i));
    end
end
