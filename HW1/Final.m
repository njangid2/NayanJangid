clear;
close all;
clc;
mu=1.327*10^11;
au=149597870.7;
day=24*60*60;

auday=au/(day);
Rearth=au*[-1.796136509111975*10^(-1) 9.667949206859814*10^(-1) -3.668681017942158*10^(-5)];
Vearth=auday*[-1.720038360888334*10^(-2) -3.211186197806460*10^(-3) 7.927736735960840*10^(-7)];

r1I = au*[3.515868886595499 * 10^-2 -3.162046390773074 4.493983111703389];
v1I = auday*[-2.317577766980901 * 10^-3 9.843360903693031 * 10^-3 -1.541856855538041 * 10^-2];

r2I = au*[7.249472033259724 14.61063037906177 14.24274452216359];
v2I = auday*[-8.241709369476881 * 10^-3 -1.156219024581502*10^-2 -1.317135977481448 * 10^-2];

% Time range for departure and arrival dates
dep_start = datenum('2017-01-01');
dep_end = datenum('2020-07-31');
arr_startb = datenum('2019-06-01');
arr_endb = datenum('2022-01-31');
arr_starto = datenum('2017-08-01');
arr_endo = datenum('2019-01-31');

% Number of points in each direction
ndep = dep_end-dep_start;
narrb = arr_endb-arr_startb;
narro = arr_endo-arr_starto;


% Create departure and arrival date vectors
dep_dates = linspace(dep_start-dep_start, dep_end-dep_start, ndep+1);
arr_datesb = linspace(arr_startb-dep_start, arr_endb-dep_start, narrb+1);
arr_dateso = linspace(arr_starto-dep_start, arr_endo-dep_start, narro+1);


for i=1:1:ndep
    [REarth(i,:) ,VEarth(i,:)] =rv_from_r0v0(Rearth, Vearth,(dep_dates(i))*day,mu);
    
end

for j=1:1:narro
    [ROum(j,:) ,VOum(j,:)] =rv_from_r0v0(r1I, v1I, (arr_dateso(j))*day,mu);
end

for j=1:1:narrb
    [Rb(j,:) ,Vb(j,:)] =rv_from_r0v0(r2I, v2I, (arr_datesb(j))*day,mu);
end



plot3(REarth(:,1),REarth(:,2),REarth(:,3))
hold on
plot3(Rb(:,1),Rb(:,2),Rb(:,3))
plot3(ROum(:,1),ROum(:,2),ROum(:,3))

legend("Earth","Borisov","Oumouamoua")
xlabel("in Km")
ylabel("in Km")
zlabel("in Km")


% orbital elem of Oumouamoua
[ho ,eo ,RAo ,inclo ,wo ,TAo ,ao]=orbit_elem_from_rv(r1I,v1I,mu)


% orbital elem of Boriisov
[hb ,eb ,RAb ,inclb ,wb ,TAb ,ab]=orbit_elem_from_rv(r2I,v2I,mu)



r1=[5644 2830 4170];
r2=[-2240 7320 4980];
mue=398600.4418;
t=1200
[v1p v2p]=lambert(r1,r2,t,'pro',mue)
[v1r v2r]=lambert(r1,r2,t,'retro',mue)