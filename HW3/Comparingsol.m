close all; 
clear all;
clc;

%...Initial orbital parameters (given):
w=[0.02 ;0.1 ;0.5];
for o=1:1:20
    if o==1
        RA0(o) = 0; %Right ascension of the node (radians)
        i0(o) = 45*pi/180; %Inclination (radians)
        w0(o) = 0; %Argument of perigee (radians)
        M0(o) = 0; %Mean anomaly (radians)
        e0(o) = 0.5; %eccentricity
        a0(o) = 1; %Semimajor axis
    else
        RA0(o) =RA0(o-1)+ 0.05; %Right ascension of the node (radians)
        i0(o) =i0(o-1)+ 0.05; %Inclination (radians)
        w0(o) =w0(o-1)+ 0.05; %Argument of perigee (radians)
        M0(o) =M0(o-1)+ 0.05; %Mean anomaly (radians)
        e0(o) =e0(o-1)+ 0.01; %eccentricity
        a0(o) =a0(o-1)+0.01; %Semimajor axis
    end
end


mu=1;
for o=1:1:20
    for u=1:1:3
        %Initial value of old variable L0 G0 H0 l0 g0 h0
        n0(o)=sqrt(1/a0(o)^3);
        L0(o)=n0(o)*a0(o)^2;
        G0(o)=L0(o)*sqrt(1-e0(o)^2);
        H0(o)=G0(o)*cos(i0(o));
        l0(o)=M0(o);
        g0(o)=w0(o);
        h0(o)=RA0(o);
        
        %Initial value of new variable L00 G00 H00 l00 g00 h00
    
        L00(o,u)=L0(o)-(-w(u)*H0(o)*L0(o)^3);%
        G00(o,u)=G0(o);%
        H00(o,u)=H0(o)-(-w(u)*L0(o)^3*l0(o));%
        l00(o,u)=l0(o)-(3*w(u)*H0(o)*L0(o)^2*l0(o));%
        g00(o,u)=g0(o);%
        h00(o,u)=h0(o)-(w(u)*L0(o)^3*l0(o));%
    end

end

for u=1:1:3
    for o=1:1:20
        %...Store initial elements (from above) in the vector coe0 :
        coe2 = [L0(o) G0(o) H0(o) l0(o) g0(o) h0(o) w(u)];
        coe1 = [l00(o,u) L00(o,u)];%no g h because they are constant
        
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
        L2=q2(:,1);
        G2=q2(:,2);
        H2=q2(:,3);
        l2 = q2(:,4);
        g2 = q2(:,5);
        h2 = q2(:,6);
        
        for ii=1:length(l2)
            hh2(ii,o,u)=e0(o)*sin(g2(ii)+h2(ii));%
            kk2(ii,o,u)=e0(o)*cos(g2(ii)+h2(ii));%
            pp2(ii,o,u)=tan(i0(o)/2)*sin(h2(ii));%
            qq2(ii,o,u)=tan(i0(o)/2)*cos(h2(ii));%
        end
        
        %Question 1 part
        ll1 = q1(:,1);
        
        
        %transforming back to old variables
        for s=1:length(ll1)
            L1(s,o,u)=L00(o,u)-(w(u)*H00(o,u)*L00(o,u)^3);
            G1(s,o,u)=G00(o,u);
            H1(s,o,u)=H00(o,u)-(w(u)*L00(o,u)^3*ll1(s));
            l1(s,o,u)=ll1(s)+(3*w(u)*H00(o,u)*L00(o,u)^2*ll1(s));
            g1(s,o,u)=g00(o,u);
            h1(s,o,u)=h00(o,u)+(w(u)*L00(o,u)^3*ll1(s));
        end
        
        
        for ii=1:length(ll1)
            
            hh1(ii,o,u)=e0(o)*sin(g1(ii,o,u)+h1(ii,o,u));
            kk1(ii,o,u)=e0(o)*cos(g1(ii,o,u)+h1(ii,o,u));
            pp1(ii,o,u)=tan(i0(o)/2)*sin(h1(ii,o,u));
            qq1(ii,o,u)=tan(i0(o)/2)*cos(h1(ii,o,u));
        end
    end
end
% change the value of u from 1 to 3 to get graph for w=0.01;0.1;0.5
u=3;
subplot(2,2,1)
plot(hh2(:,:,u),kk2(:,:,u))
xlabel("h")
ylabel("k")
title("Plot from Hamiltonian expressed in Delaunay variables (Question2)")

subplot(2,2,2)
plot(hh1(:,:,u),kk1(:,:,u))
xlabel("h")
ylabel("k")
title("Plot from Hamiltonian expressed in New Delaunay variables (Question1)")

subplot(2,2,3)
plot(pp2(:,:,u),qq2(:,:,u))
xlabel("p")
ylabel("q")
title("Plot from Hamiltonian expressed in Delaunay variables (Question2)")

subplot(2,2,4)
plot(pp1(:,:,u),qq1(:,:,u))
xlabel("p")
ylabel("q")
title("Plot from Hamiltonian expressed in New Delaunay variables (Question1)")





function dfdt = rates2(t,f)
% equation of motions (eq 31-33 from report) 
L=f(1);
G=f(2);
H=f(3);
l = f(4);
g = f(5);
h = f(6);
w=f(7);%rotation rate
Ldot=0;
Gdot=0;
Hdot=0;
ldot=1/(L^3);
gdot=0;
hdot=w;
wdot=0;

dfdt = [Ldot Gdot Hdot ldot gdot hdot wdot]';
end 

function dfdt = rates1(t,f)
% equation of motions (eq 31-33 from report) 
l = f(1);
L=f(2);
ldot=1/(L^3);
Ldot=0;


dfdt = [ldot Ldot]';
end 








