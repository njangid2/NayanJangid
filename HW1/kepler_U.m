function x = kepler_U(dt, ro, vro, a,mu)
    error =1.e-8;
    nMax =1000;
    x=sqrt(mu)*abs(a)*dt;
    n=0;
    ratio =1;
    while abs(ratio) > error && n <=nMax
        n=n+1;
        C=stumpC(a*x^2);
        S=stumpS(a*x^2);
        F=(ro*vro/sqrt(mu))*x^2*C + (1 - a*ro)*x^3*S + ro*x - sqrt(mu)*dt;
        dFdx =(ro*vro/sqrt(mu))*x*(1 - a*x^2*S) + (1 - a*ro)*x^2*C + ro;
        ratio =F/dFdx;
        x=x - ratio;
    end
end