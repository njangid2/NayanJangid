clear;
close all;
clc;
mu=1.327*10^11;
au=149597870.7;
day=24*60*60;

auday=au/(day);
Rearth=au*[-1.796136509111975*10^(-1) 9.667949206859814*10^(-1) -3.668681017942158*10^(-5)];
Vearth=auday*[-1.720038360888334*10^(-2) -3.211186197806460*10^(-3) 7.927736735960840*10^(-7)];

r2I = au*[7.249472033259724 14.61063037906177 14.24274452216359];
v2I = auday*[-8.241709369476881 * 10^-3 -1.156219024581502*10^-2 -1.317135977481448 * 10^-2];

% Time range for departure and arrival dates
dep_start = datenum('2017-01-01');
dep_end = datenum('2020-07-31');
arr_start = datenum('2019-06-01');
arr_end = datenum('2022-01-31');

% Number of points in each direction
ndep = dep_end-dep_start;
narr = arr_end-arr_start;


% Create departure and arrival date vectors
dep_dates = linspace(dep_start-dep_start, dep_end-dep_start, ndep+1);
arr_dates = linspace(arr_start-dep_start, arr_end-dep_start, narr+1);

dv_rendezvous = zeros(narr, ndep);
dv_flyby = zeros(narr, ndep);
 
for i=1:1:ndep
    [REarth(i,:) ,VEarth(i,:)] =rv_from_r0v0(Rearth, Vearth,(dep_dates(i))*day,mu);
    
end
for j=1:1:narr
    [RB(j,:) ,VB(j,:)] =rv_from_r0v0(r2I, v2I, (arr_dates(j))*day,mu);
end

for i=1:ndep
    re=REarth(i,:);
    ve=VEarth(i,:);
    for j=1:narr
          if arr_dates(j)>dep_dates(i)
            dt = (arr_dates(j) - dep_dates(i)) * day;
            rb=RB(j,:);
            vb =VB(j,:);
            [v1, v2]=lambert(re,rb,dt,'pro',mu); 
            a=norm(v1 - ve) + norm(vb - v2);
            b=norm(v1 - ve);
            if a<=60
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
xticks([0 59 120 181 243 304 365 424 485 546 608 669 730 789 850 911 973 1034 1095 1155 1216 1277])
xticklabels(["1-Jan-17","1-Mar-17","1-May-17","1-Jul-17","1-Sep-17","1-Nov-17","1-Jan-18","1-Mar-18","1-May-18","1-Jul-18","1-Sep-18","1-Nov-18","1-Jan-19","1-Mar-19","1-May-19","1-Jul-19","1-Sep-19","1-Nov-19","1-Jan-20","1-Mar-20","1-May-20","1-Jul-20"])
ylabel('Arrival date')
yticks([1 61 122 183 245 305 366 427 488 549 611 670 731 792 853 914])
yticklabels(["1-Jun-19","1-Aug-19","1-Oct-19","1-Dec-19","1-Feb-20","1-Apr-20","1-Jun-20","1-Aug-20","1-Oct-20","1-Dec-20","1-Feb-21","1-Apr-21","1-Jun-21","1-Aug-21","1-Oct-21","1-Dec-21"])
subplot(2,1,2)
contourf(dv_flyby)   
colorbar
title("fly-by mission")
xlabel('Departure date')
xticks([0 59 120 181 243 304 365 424 485 546 608 669 730 789 850 911 973 1034 1095 1155 1216 1277])
xticklabels(["1-Jan-17","1-Mar-17","1-May-17","1-Jul-17","1-Sep-17","1-Nov-17","1-Jan-18","1-Mar-18","1-May-18","1-Jul-18","1-Sep-18","1-Nov-18","1-Jan-19","1-Mar-19","1-May-19","1-Jul-19","1-Sep-19","1-Nov-19","1-Jan-20","1-Mar-20","1-May-20","1-Jul-20"])
ylabel('Arrival date')
yticks([1 61 122 183 245 305 366 427 488 549 611 670 731 792 853 914])
yticklabels(["1-Jun-19","1-Aug-19","1-Oct-19","1-Dec-19","1-Feb-20","1-Apr-20","1-Jun-20","1-Aug-20","1-Oct-20","1-Dec-20","1-Feb-21","1-Apr-21","1-Jun-21","1-Aug-21","1-Oct-21","1-Dec-21"])