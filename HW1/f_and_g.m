function [f, g] =f_and_g(x, t, ro, a,mu)
    z=a*x^2;
    f=1 - (x^2/ro)*stumpC(z);
    g=t - (1/sqrt(mu))*x^3*stumpS(z);
end