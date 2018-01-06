function approx = APPR_tB_baryc_two(x, lambda, sample, func, weight, Cd)
    a  = sample;
    n  = length(a);
    f  = func;
    w  = weight;
    Sd = @(z) (1/2)*sinh(2*Cd*z);
    
    for i=1:n
        R(i,:) = lambda(i)./Sd(x - a(i));
    end

    fwsum = 0;
    for i=1:n
        fwsum = fwsum + R(i,:) * f(a(i))./w(a(i));
    end
    
    approx = w(x) .* fwsum ./ sum(R);
end