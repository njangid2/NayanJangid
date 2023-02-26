function [R,V] =rv_from_r0v0(R0, V0, t,mu)
    r0 =norm(R0);
    v0 =norm(V0);
    vr0 =dot(R0, V0)/r0;
    alpha =2/r0 - v0^2/mu;
    x=kepler_U(t, r0, vr0, alpha,mu);
    [f, g] =f_and_g(x, t, r0, alpha,mu);
    R=f*R0 + g*V0;
    r=norm(R);
    [fdot, gdot] =fDot_and_gDot(x, r, r0, alpha,mu);
    V=fdot*R0 + gdot*V0;
end