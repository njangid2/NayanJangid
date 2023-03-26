close all; 
clear all;
clc;

% Q1

%INITIAL VALUES

Te=8*60*60; %time period
J2=0.00108;
R=6370;
mu=398600;
a=(((Te/(2*pi))^2)*mu)^(1/3);
n=sqrt(mu/(a^3));

i=[acos(sqrt(1/5)) acos(-sqrt(1/5))] %solution in range of [0,pi]
j=1;
for r=600:50:a
    e(j)=1-((R+r)/a);
    
    for u=1:2
        DR(j,u)=-(3/2)*n*J2*((R/a)^2)*((cos(i(u)))/((1-(e(j)^2))^2));
    end
    j=j+1;
end
r=600:50:a;
plot(r,DR(:,1))
hold on
plot(r,DR(:,2))
hold off
legend("i = 1.1071","i = 2.0344")
xlabel("Altitude (Km)")
ylabel("Drift rate")
title("Drift rate vs Perigee Altitude")
u=1;
j=1;
% for u=1 and j=1 DR is lowest
lowest_driftrate=(DR(j,u))
semimajoraxis=vpa(a)
incli=vpa(i(u))
ecc=vpa(e(j))


