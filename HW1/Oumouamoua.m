clear all;
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

% Time range for departure and arrival dates
dep_start = datenum('2017-01-01');
dep_end = datenum('2017-12-31');
arr_start = datenum('2017-08-01');
arr_end = datenum('2019-01-31');

% Number of points in each direction
ndep = dep_end-dep_start;
narr = arr_end-arr_start;


% Create departure and arrival date vectors
dep_dates = linspace(dep_start-dep_start, dep_end-dep_start, ndep+1);
arr_dates = linspace(arr_start-dep_start, arr_end-dep_start, narr+1);


% Create arrays for Rearth, Vearth, r1I, and v1I


for i=1:1:ndep
    [REarth(i,:) ,VEarth(i,:)] =rv_from_r0v0(Rearth, Vearth,(dep_dates(i))*day,mu);
    
end
for j=1:1:narr
    [ROum(j,:) ,VOum(j,:)] =rv_from_r0v0(r1I, v1I, (arr_dates(j))*day,mu);
end
dv_rendezvous = zeros(narr, ndep);
dv_flyby = zeros(narr, ndep);
for i=1:ndep
    re=REarth(i,:);
    ve=VEarth(i,:);
    for j=1:narr
          if arr_dates(j)>dep_dates(i)
            dt = (arr_dates(j) - dep_dates(i)) * day;
            rb=ROum(j,:);
            vb =VOum(j,:);
            [v1, v2]=lambert(re,rb,dt,'pro',mu); 
            a=norm(v1 - ve) + norm(vb - v2);
            b=norm(v1 - ve);
            if a<=50
            dv_rendezvous(j, i) =norm(v1 - ve) + norm(vb - v2) ;
            end
            if b<=20
            dv_flyby(j, i) =norm(v1 - ve) ;
            end
          end
    end
end
dv_rendezvous(dv_rendezvous == 0) = NaN;
dv_flyby(dv_flyby == 0) = NaN;
subplot(2,1,1)
contourf(dv_rendezvous)
colorbar
title("Rendezvous mission")
xlabel('Departure date')
xticks([0 31 59 90 120 151 181 212 243 273 304 334])
xticklabels(["1-Jan-17","1-Feb-17","1-Mar-17","1-Apr-17","1-May-17","1-Jun-17","1-Jul-17","1-Aug-17","1-Sep-17","1-Oct-17","1-Nov-17","1-Dec-17"])
ylabel('Arrival date')
yticks([1 31 61 92 122 153 184 213 243 273 304 334 365 396 426 457 487 518])
yticklabels(["1-Aug-17","1-Sep-17","1-Oct-17","1-Nov-17","1-Dec-17","1-Jan-18","1-Feb-18","1-Mar-18","1-Apr-18","1-May-18","1-Jun-18","1-Jul-18","1-Aug-18","1-Sep-18","1-Oct-18","1-Nov-18","1-Dec-18","1-Jan-19"])
subplot(2,1,2)
contourf(dv_flyby)   
colorbar
title("fly-by mission")
xlabel('Departure date')
xticks([0 31 59 90 120 151 181 212 243 273 304 334])
xticklabels(["1-Jan-17","1-Feb-17","1-Mar-17","1-Apr-17","1-May-17","1-Jun-17","1-Jul-17","1-Aug-17","1-Sep-17","1-Oct-17","1-Nov-17","1-Dec-17"])
ylabel('Arrival date')
yticks([1 31 61 92 122 153 184 213 243 273 304 334 365 396 426 457 487 518])
yticklabels(["1-Aug-17","1-Sep-17","1-Oct-17","1-Nov-17","1-Dec-17","1-Jan-18","1-Feb-18","1-Mar-18","1-Apr-18","1-May-18","1-Jun-18","1-Jul-18","1-Aug-18","1-Sep-18","1-Oct-18","1-Nov-18","1-Dec-18","1-Jan-19"])