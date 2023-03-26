close all; 
clear all;
clc;
%...Conversion factors:
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180; %Degrees to radians

%...Constants:
mu = 398600; %Gravitational parameter (km^3/s^2)
RE = 6378; %Earth’s radius (km)
J2 = 0.00108; %Earth*s J2

%...Initial orbital parameters (given):
H0 = 0; 
i0 = 0; 
p0 = 0; 
q0 = 0; 
lam0 = 0;
e0 = 0.05; 
a0 = 7000; 
k0 = e0;
M0=0;
TA0=kepler_TA(e0,M0);
L0=TA0; % w + raan =0

T0 = 2*(pi/sqrt(mu))*a0^1.5; %Period (s)
%...Store initial orbital elements (from above) in the vector coe0:
coe0 = [a0 H0 k0 p0 q0 L0];

%...Use ODE45 to integrate the Gauss variational equations (Equations
% 12.89) from t0 to tf:
t0 = 0;
tf = 2*days;
nout = 5000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
options = odeset('reltol', 1.e-8, 'abstol', 1.e-8, 'initialstep', T0/1000);
y0 = coe0';
[t,y] = ode45(@rates, tspan, y0, options);
%Assign the time histories mnemonic variable names:
a = y(:,1);
H = y(:,2);
k = y(:,3);
p = y(:,4);
q = y(:,5);
L = y(:,6);

%Plot 
figure(1)
subplot(5,1,1)
plot(t/days,a)
xlabel('days')
grid on
grid minor
axis tight
title("a in km")

subplot(5,1,2)
plot(t/days,H)
xlabel('days')
grid on
grid minor
axis tight
title("h")

subplot(5,1,3)
plot(t/days,k)
xlabel('days')
grid on
grid minor
axis tight
title("k")

subplot(5,1,4)
plot(t/days,p )
xlabel('days')
grid on
grid minor
axis tight
title("p")

subplot(5,1,5)
plot(t/days,q)
xlabel('days')
grid on
grid minor
axis tight
title("q")





function dfdt = rates(t,f)
% This function calculates the time rates of the orbital elements
%...The orbital elements at time t:

mu = 398600; %Gravitational parameter (km^3/s^2)
RE = 6378; %Earth’s radius (km)
J2 = 1082.63e-6; %Earth’s J2
a = f(1);
h = f(2);
k = f(3);
p = f(4);
q = f(5);
L = f(6);
B=sqrt(1-h^2-k^2);
phi=1+h*sin(L)+k*cos(L);
G=1+p^2+q^2;
ft=((12*mu*J2*RE^2)/((B^8)*(G^2)*(a^4)))*(q*cos(L)+p*sin(L))*(p*cos(L)-q*sin(L))*phi^4;
fh=((6*mu*J2*RE^2)/((B^8)*(G^2)*(a^4)))*(p*cos(L)-q*sin(L))*(1-p^2-q^2)*phi^4;
fr=((3*mu*J2*RE^2)/(2*(B^8)*(a^4)))*((12/(G^2))*((p*cos(L)-q*sin(L))^2)-1)*phi^4;
% Orbital element rates at time t (Equations 12.89)(curtis):
adot=(2/B)*sqrt((a^3)/mu)*((k*sin(L)-h*cos(L))*fr+phi*ft);
hdot=B*sqrt(a/mu)*(-fr*cos(L)+(((h+sin(L))/phi)+sin(L))*ft-(k/phi)*(p*cos(L)-q*sin(L))*fh);
kdot=B*sqrt(a/mu)*(fr*sin(L)+(((k+cos(L))/phi)+cos(L))*ft+(h/phi)*(p*cos(L)-q*sin(L))*fh);
pdot=(B/2)*sqrt(a/mu)*(1+p^2+q^2)*(sin(L)/phi)*fh;
qdot=(B/2)*sqrt(a/mu)*(1+p^2+q^2)*(cos(L)/phi)*fh;
Ldot=((phi^2)/(B^3))*sqrt(mu/(a^3));
dfdt = [adot hdot kdot pdot qdot Ldot]';
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



