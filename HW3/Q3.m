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

%Initial value of old variable L0 G0 H0 l0 g0 h0
L0=n0*a0^2;

% G0=-0.8660
G0=L0*sqrt(1-e0^2);
H0=G0*cos(i0);
l0=M0;
g0=w0;
h0=RA0;
%rotation rate
w=sqrt(0.02^2+0.1^2+0.5^2);

%Initial value of new variable L00 G00 H00 l00 g00 h00
L00=L0-(-w*H0*L0^3);
G00=G0;
H00=H0-(-w*L0^3*l0);
l00=l0-(3*w*H0*L0^2*l0);
g00=g0;
h00=h0-(w*L0^3*l0);



%...Store initial elements (from above) in the vector coe0 ( no L, G, H because they are constant and hence a0, e0 and i0 is constant):
coe2 = [l0 g0 h0];
coe1 = [l00 L00];%no g h because they are constant

%...Use ODE45 to integrate from t0 to tf:
t0 = 0;
tf = 100;
nout = 2000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
options = odeset('reltol', 1.e-8, 'abstol', 1.e-8);
y2 = coe2';
y1 = coe1';
[t,q2] = ode45(@rates2, tspan, y2, options);
[t,q1] = ode45(@rates1, tspan, y1, options);

%Question2 part
%Assign the time histories mnemonic variable names:
l2 = q2(:,1);
g2 = q2(:,2);
h2 = q2(:,3);

for ii=1:length(l2)
    hh2(ii)=e0*sin(g2(ii)+h2(ii));
    kk2(ii)=e0*cos(g2(ii)+h2(ii));
    pp2(ii)=tan(i0/2)*sin(h2(ii));
    qq2(ii)=tan(i0/2)*cos(h2(ii));
end

%Question 1 part
ll1 = q1(:,1);


%transforming back to old variables
for o=1:length(ll1)
    L1(o)=L00-(w*H00*L00^3);
    G1(o)=G00;
    H1(o)=H00-(w*L00^3*ll1(o));
    l1(o)=ll1(o)+(3*w*H00*L00^2*ll1(o));
    g1(o)=g00;
    h1(o)=h00+(w*L00^3*ll1(o));
end


for ii=1:length(ll1)
    
    hh1(ii)=e0*sin(g1(ii)+h1(ii));
    kk1(ii)=e0*cos(g1(ii)+h1(ii));
    pp1(ii)=tan(i0/2)*sin(h1(ii));
    qq1(ii)=tan(i0/2)*cos(h1(ii));
end

subplot(2,2,1)
plot(hh2,kk2)
xlabel("h")
ylabel("k")
title("Plot from Hamiltonian expressed in Delaunay variables (Question2)")

subplot(2,2,2)
plot(hh1,kk1)
xlabel("h")
ylabel("k")
title("Plot from Hamiltonian expressed in New Delaunay variables (Question1)")

subplot(2,2,3)
plot(pp2,qq2)
xlabel("p")
ylabel("q")
title("Plot from Hamiltonian expressed in Delaunay variables (Question2)")

subplot(2,2,4)
plot(pp1,qq1)
xlabel("p")
ylabel("q")
title("Plot from Hamiltonian expressed in New Delaunay variables (Question1)")





function dfdt = rates2(t,f)
% equation of motions (eq 31-33 from report) 
l = f(1);
g = f(2);
h = f(3);
w=sqrt(0.02^2+0.1^2+0.5^2);%rotation rate
L=1;%L0
i = 45*pi/180;
e = 0.5;
ldot=1/(L^3);
gdot=0;
hdot=w;

dfdt = [ldot gdot hdot]';
end 

function dfdt = rates1(t,f)
% equation of motions (eq 31-33 from report) 
l = f(1);
L=f(2);
ldot=1/(L^3);
Ldot=0;


dfdt = [ldot Ldot]';
end 




