close all; 
clear all;
clc;

%...Initial orbital parameters (given):
RA0 = 0; %Right ascension of the node (radians)
i0 = 45*pi/180; %Inclination (radians)
w0 = 0; %Argument of perigee (radians)
M0 = 0; %Mean anomaly (radians)
e0 = 0.5; %eccentricity
a0 = 1; %Semimajor axis
n0=sqrt(1/a0^3);
L0=n0*a0^2;
mu=1;
% G0=-0.8660
G0=L0*sqrt(1-e0^2);
H0=G0*cos(i0);
l0=M0;
g0=w0;
h0=RA0;


%...Store initial elements (from above) in the vector coe0 ( no L, G, H because they are constant and hence a0, e0 and i0 is constant):
coe0 = [l0 g0 h0];
angular=sqrt(mu*a0*(1-e0)^2);
%...Use ODE45 to integrate from t0 to tf:
t0 = 0;
tf = 100;
nout = 2000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
options = odeset('reltol', 1.e-8, 'abstol', 1.e-8);
y0 = coe0';
[t,y] = ode45(@rates, tspan, y0, options);
%Assign the time histories mnemonic variable names:
l = y(:,1);
g = y(:,2);
h = y(:,3);
length(l);
for ii=1:length(l)
    TA(ii)=kepler_TA(e0,l(ii));
    [r(ii,:),v(ii,:)]=sv_from_coe([angular e0 h(ii) i0 g(ii) TA(ii)],mu);
end

plot3(r(:,1),r(:,2),r(:,3))
hold on
plot3(0,0,0,"*")
hold off
xlabel("x in DU")
ylabel("y in DU")
zlabel("z in DU")
title("Orbit")






function dfdt = rates(t,f)
% equation of motions (eq 31-33 from report) 
l = f(1);
g = f(2);
h = f(3);
w=0.01;%rotation rate
L=1;%L0
i = 45*pi/180;
e = 0.5;
ldot=1/(L^3);
gdot=0;
hdot=w;

dfdt = [ldot gdot hdot]';
end 

function f = kepler_TA(e, M)

%given the eccentricity and the mean anomaly.
%E - eccentric anomaly (radians)
%e - eccentricity, passed from the calling program
%M - mean anomaly (radians), passed from the calling program
%f - true anomaly
%.Select a starting value for E:
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end
%Iterate on Equation 3.17(curtis) until E is determined to within
%the error tolerance:
ratio = 1;
while abs(ratio) > 10^(-10)
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end
f=2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));
end 



function [r, v] = sv_from_coe(coe,mu)
h = coe(1);% angular momentum
e = coe(2);% ecc
RA = coe(3);% RAAN
incl = coe(4);%incl
w = coe(5);%argument of periapsis
TA = coe(6);% true ano
%...Equations 4.45 and 4.46 (rp and vp are column vectors) Curtis
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);
%...Equation 4.34:Curtis
R3_W = [ cos(RA) sin(RA) 0;-sin(RA) cos(RA) 0;0 0 1];

%...Equation 4.32:Curtis
R1_i = [1 0 0;0 cos(incl) sin(incl);0 -sin(incl) cos(incl)];
%...Equation 4.34:Curtis
R3_w = [ cos(w) sin(w) 0;-sin(w) cos(w) 0;0 0 1];
%...Equation 4.49:Curtis
Q_pX = (R3_w*R1_i*R3_W)';
%...Equations 4.51 (r and v are column vectors):Curtis
r = Q_pX*rp;
v = Q_pX*vp;
%...Convert r and v into row vectors:
r = r';
v = v';
end


