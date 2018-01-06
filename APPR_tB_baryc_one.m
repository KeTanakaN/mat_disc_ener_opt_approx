function approx = APPR_tB_baryc_one(x, lambda, sample, func, weight, Cd)
    a  = sample;
    n  = length(a);
    f  = func;
    w  = weight;
    Sd = @(z) (1/2)*sinh(2*Cd*z);
    
    tB = 1;
    for i=1:n
        tB = tB .* tanh(Cd*(x - a(i)));
    end

    fwsum = 0;
    for i=1:n
        fwsum = fwsum + lambda(i)*f(a(i))./(Sd(x - a(i))*w(a(i)));
    end
    
    approx = w(x) .* tB .* fwsum;
end