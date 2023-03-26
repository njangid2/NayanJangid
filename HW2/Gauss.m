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
RA0 = 90*deg; %Right ascension of the node (radians)
i0 = 1.10654; %Inclination (radians)
w0 = 5*deg; %Argument of perigee (radians)
M0 = 10*deg; %Mean anomaly (radians)
e0 = 0.74; %eccentricity
a0 = 26600; %Semimajor axis (km)

h0 = sqrt(a0*mu*(1 - e0^2)); %angular momentum (km^2/s)
T0 = 2*(pi/sqrt(mu))*a0^1.5; %Period (s)
TA0=kepler_TA(e0,M0); %True anomaly (radians)
%...Store initial orbital elements (from above) in the vector coe0:
coe0 = [h0 e0 RA0 i0 w0 TA0];

%...Use ODE45 to integrate the Gauss variational equations (Equations
% 12.89) from t0 to tf:
t0 = 0;
tf = 100*days;
nout = 5000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
options = odeset('reltol', 1.e-8, 'abstol', 1.e-8, 'initialstep', T0/1000);
y0 = coe0';
[t,y] = ode45(@rates, tspan, y0, options);
%Assign the time histories mnemonic variable names:
h = y(:,1);
e = y(:,2);
RA = y(:,3);
i = y(:,4);
w = y(:,5);
TA = y(:,6);
for ii=1:length(h)
    a(ii)=(h(ii)^2)/(mu*(1-e(ii)^2));
end
%Plot 
figure(1)
subplot(5,1,1)
plot(t/days,RA /deg)
title('Right Ascension (degrees)')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,2)
plot(t/days,w /deg)
title('Argument of Perigee (degrees)')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,3)
plot(t/days,a)
title('semimajor axis (km)')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,4)
plot(t/days,e )
title('Eccentricity')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,5)
plot(t/days,i /deg)
title('Inclination (degrees)')
xlabel('days')
grid on
grid minor
axis tight




function dfdt = rates(t,f)
% This function calculates the time rates of the orbital elements using Gaussss variational equations (Equations 12.89)(curtis).
%...The orbital elements at time t:

mu = 398600; %Gravitational parameter (km^3/s^2)
RE = 6378; %Earth’s radius (km)
J2 = 1082.63e-6; %Earth’s J2
h = f(1);
e = f(2);
RA = f(3);
i = f(4);
w = f(5);
TA = f(6);
r = h^2/mu/(1 + e*cos(TA)); %The radius
u = w + TA; %Argument of latitude
% Orbital element rates at time t (Equations 12.89)(curtis):
hdot = -(3/2)*((J2*mu*RE^2)/(r^3))*((sin(i))^2)*sin(2*u);
edot = (3/2)*((J2*mu*RE^2)/(h*r^3))*(((h^2)/(mu*r))*sin(TA)*(3*sin(i)^2*sin(u)^2 - 1)-sin(2*u)*sin(i)^2*((2+e*cos(TA))*cos(TA)+e));
TAdot = h/r^2 + (3/2)*((J2*mu*RE^2)/(e*h*r^3))*(((h^2)/(mu*r))*cos(TA)*(3*sin(i)^2*sin(u)^2 - 1)+ sin(2*u)*sin(i)^2*sin(TA)*(e*cos(TA) + 2));
RAdot = -3*((J2*mu*RE^2)/(h*r^3))*(sin(u)^2)*cos(i);
idot = (-3/4)*((J2*mu*RE^2)/(h*r^3))*sin(2*u)*sin(2*i);
wdot = (3/2)*((J2*mu*RE^2)/(e*h*r^3))*(((h^2)/(mu*r))*cos(TA)*(1-3*sin(i)^2*sin(u)^2)- sin(2*u)*sin(i)^2*sin(TA)*(2 + e*cos(TA))+ 2*e*cos(i)^2*sin(u)^2);
dfdt = [hdot edot RAdot idot wdot TAdot]';
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


