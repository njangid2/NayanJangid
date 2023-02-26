function [fdot, gdot] =fDot_and_gDot(x,r,ro,a,mu)
    z=a*x^2;
    fdot =(sqrt(mu)/(r*ro))*(z*stumpS(z) - 1)*x;
    gdot =1 - (x^2/r)*stumpC(z);
end