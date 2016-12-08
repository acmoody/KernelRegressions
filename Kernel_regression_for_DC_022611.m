clear all
close all

% sitecode key
afnames(1,:) = 'US-Seg'; % 1-GLand
afnames(2,:) = 'US-Ses'; % 2-SLand
afnames(3,:) = 'US-Wjs'; % 3-JSav
afnames(4,:)='US-Mpj'; % 4-PJ
afnames(5,:)='US-Vcp'; % 5-PPine
afnames(6,:)='US-Vcm'; % 6-MCon
afnames(7,:)='US-FR2'; % 7-TX_savanna

colour(1,:)=[0.9 0.5 0.0];
colour(2,:)=[0.6 0.2 0];
colour(3,:)=[0.25 1.0 0.0];
colour(4,:)=[0.0 0.5 0.0];
colour(5,:)=[0.5 0.5 1.0];
colour(6,:)=[0.0 0.0 0.6];

firstday(1)=151;
firstday(2)=90;
firstday(3)=59;
firstday(4)=59;
firstday(5)=90;
firstday(6)=120;

lastday(1)=272;
lastday(2)=272;
lastday(3)=303;
lastday(4)=303;
lastday(5)=303;
lastday(6)=303;

figure;
shallow=gcf;
figure;
deep=gcf;
figure;
deep_n=gcf;
figure;
shallow_n=gcf;


for jj = 1:6

%%
sitecode=jj;

%%
if sitecode==1
year = 2007;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2);
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data1=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2008;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data2=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2010;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/US-Sen_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data4=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

data=cat(1,data1,data2,data4);

elseif sitecode ==2
year = 2007;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2);
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data1=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2008;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data2=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2009;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data3=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

data=cat(1,data1,data2,data3);

elseif sitecode == 3

year = 2007;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2);
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data1=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2008;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data2=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2009;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,47);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data3=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2010;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,47);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data4=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

data=cat(1,data1,data2,data3,data4);

elseif sitecode == 4
%%
year = 2008;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data2=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2009;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data3=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2010;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data4=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

data=cat(1,data2,data3,data4);

elseif sitecode==5 || sitecode ==6

year = 2007;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2);
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data1=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2008;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data2=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2009;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data3=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

year = 2010;
year_s=num2str(year);
filename = strcat('../Ameriflux_files/',afnames(sitecode,:),'_',year_s,'_gapfilled.txt');

ffi=dlmread(filename,'',5,0);
ffi(ffi==-9999)=nan;

yr=ffi(:,1);
day=ffi(:,2)+365+366+365;
hr=ffi(:,3);
dt=ffi(:,4);
ta=ffi(:,6);
ppt=ffi(:,22);
swc=ffi(:,28);
nee=ffi(:,11);
swc2=ffi(:,46);
re=ffi(:,41);
gpp=ffi(:,43);
rg=ffi(:,33);
par=ffi(:,30);   
par(isnan(par))=rg(isnan(par)).*2.1032-8.2985;

data4=cat(2,yr,day,hr,dt,ta,ppt,swc,swc2,nee,gpp,re,par);

data=cat(1,data1,data2,data3,data4);
end

%%
data(data==-9999)=nan;
mu2e=(1*60*30)./1000000;
mu2g=((1./1000000)*12)*60*30;

ndays=(length(data)/48);

for i = 1:ndays
    dayy=find(data(:,2)==i+(data(1,2)-1));
    touse=data(dayy,:);
    
    out(i,1)=nanmean(touse(:,1));       % 1 year
    out(i,2)=nanmean(touse(:,2));       % 2 day 
    out(i,3)=nanmean(touse(:,3));       % 3 hour
    out(i,4)=nanmean(touse(:,4));       % 4 decitime
    out(i,5)=nanmean(touse(:,5));       % 5 air temp
    out(i,6)=sum(touse(:,6));           % 6 precip
    out(i,7)=nanmean(touse(:,7));       % 7 swc shallow
    out(i,8)=nanmean(touse(:,8));       % 8 swc deep
    out(i,9)=(sum(touse(:,9))).*mu2g;   % 9 nee
    out(i,10)=(sum(touse(:,10))).*mu2g; % 10 gpp
    out(i,11)=(sum(touse(:,11))).*mu2g; % 11 re
    
    daytime=find(touse(:,3)>530 & touse(:,3)<1830);
    nigtime=find(touse(:,3)<600 | touse(:,3)>1800);
    
    out(i,12)=(sum(touse(daytime,9))).*mu2g;   % 12 daytime nee
    out(i,13)=(sum(touse(nigtime,9))).*mu2g;   % 13 nighttime nee
    out(i,14)=nanmean(touse(daytime,12));   % Mean of daytime PAR
    if i == 1
        out(i,15)=0;
    elseif out(i,6)>0
        out(i,15)=0;
    else
        out(i,15)=out(i-1,15)+1;
    end
    
end


%% Temperature response

% dayrange=find(firstday(sitecode) < out(:,4) & out(:,4)< lastday(sitecode));
dayrange=find(0 < out(:,4) & out(:,4)< 100000);

d = out(dayrange,2);
x = out(dayrange,7);
x(x<0.001)=0.001;
y = out(dayrange,5);
z = out(dayrange,13);

x = log(x);
lt=-10; % minimum temp
tr=40; % temp range
ld=log(0.01);
dr=log(0.38)-log(0.01);
tsteps=10; % number of steps
dsteps=10;
cnt=0;

for i = 1:tsteps
    mint=(lt)+((i-1).*(tr./tsteps));
    maxt=(lt)+i.*(tr./tsteps);
    valut=(mint+maxt)./2;
    for j = 1:dsteps
        cnt=cnt+1;
        mind=(ld)+((j-1).*(dr./dsteps));
        maxd=(ld)+j.*(dr./dsteps);
        valud=(mind+maxd)./2;    
        XX2(i)=valut; YY2(j)=valud;
        found=find((y>=mint & y<maxt)&(x>=mind & x<maxd));
        flux(j,i)=nanmean(z(found));
        n_count(j,i)=0;
        n_count(j,i)=length(found);
    end
end

[XX,YY]=meshgrid(XX2(:),YY2(:));
r=ksrmv([XX(:) YY(:)],flux(:));
fluxr=flux;
fluxr(:)=r.f;


%%
left=[0.1 0.55 0.1 0.55 0.1 0.55];
bottom=[0.7 0.7 0.4 0.4 0.1 0.1];

figure(shallow);
subplot('Position',[left(sitecode) bottom(sitecode) 0.4 0.25])
 
contourf(XX,exp(YY),fluxr); 

if sitecode==3
    ylabel('Shallow SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
end
if sitecode==5 || sitecode==6
    xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
end
set(gca,'fontweight','bold','fontsize',12)
colorbar
ymin=max(0.02,min(exp(x))); ymax=min(0.31,max(exp(x)));
xmin=max(-8,min(y)); xmax=min(28,max(y));
ylim([ymin ymax]); xlim([xmin xmax]);
%ylim([0.02 0.31]); xlim([-8 28]);

n_count=n_count./(sum(reshape(n_count,100,1)));

figure(shallow_n);
subplot('Position',[left(sitecode) bottom(sitecode) 0.4 0.25])

contourf(XX,exp(YY),n_count.*fluxr); 

if sitecode==3
    ylabel('Shallow SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
end
if sitecode==5 || sitecode==6
    xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
end
set(gca,'fontweight','bold','fontsize',12)
colorbar
ymin=max(0.02,min(exp(x))); ymax=min(0.31,max(exp(x)));
xmin=max(-8,min(y)); xmax=min(28,max(y));
ylim([ymin ymax]); xlim([xmin xmax]);
%ylim([0.02 0.31]); xlim([-8 28]);

clear x; clear y; clear z;
%%

x = out(dayrange,8);
x(x<0.001)=0.001;
y = out(dayrange,5);
z = out(dayrange,12);

x = log(x);
lt=-10; % minimum temp
tr=40; % temp range
ld=log(0.07);
dr=log(0.38)-log(0.07);
tsteps=10; % number of steps
dsteps=10;
cnt=0;

for i = 1:tsteps
    mint=(lt)+((i-1).*(tr./tsteps));
    maxt=(lt)+i.*(tr./tsteps);
    valut=(mint+maxt)./2;
    for j = 1:dsteps
        cnt=cnt+1;
        mind=(ld)+((j-1).*(dr./dsteps));
        maxd=(ld)+j.*(dr./dsteps);
        valud=(mind+maxd)./2;    
        XX2(i)=valut; YY2(j)=valud;
        found=find((y>=mint & y<maxt)&(x>=mind & x<maxd));
        flux(j,i)=nanmean(z(found));
        n_count(j,i)=0;
        n_count(j,i)=length(found);
    end
end

[XX,YY]=meshgrid(XX2(:),YY2(:));
r=ksrmv([XX(:) YY(:)],flux(:));
fluxr=flux;
fluxr(:)=r.f;

% Subsample grid
XX3=linspace(XX2(1),XX2(length(XX2)),10);
YY3=linspace(YY2(1),YY2(length(YY2)),10);
[XX4,YY4]=meshgrid(XX3(:),YY3(:));
ZI = interp2(XX,YY,fluxr,XX4,YY4);


figure(deep);
subplot('Position',[left(sitecode) bottom(sitecode) 0.4 0.25])
 
contourf(XX,exp(YY),fluxr); 
 
if sitecode==3
    ylabel('Deep SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
end
if sitecode==5 || sitecode==6
    xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
end
set(gca,'fontweight','bold','fontsize',12)
colorbar
ymin=max(0.075,min(exp(x))); ymax=min(0.35,max(exp(x)));
xmin=max(-8,min(y)); xmax=min(28,max(y));
ylim([ymin ymax]); xlim([xmin xmax]);
%ylim([0.075 0.35]); xlim([-8 28]);


n_count=n_count./(sum(reshape(n_count,100,1)));
figure(deep_n);
subplot('Position',[left(sitecode) bottom(sitecode) 0.4 0.25])

contourf(XX,exp(YY),n_count.*fluxr); 
 
if sitecode==3
    ylabel('Deep SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
end
if sitecode==5 || sitecode==6
    xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
end
set(gca,'fontweight','bold','fontsize',12)
colorbar
ymin=max(0.075,min(exp(x))); ymax=min(0.35,max(exp(x)));
xmin=max(-8,min(y)); xmax=min(28,max(y));
ylim([ymin ymax]); xlim([xmin xmax]);
%ylim([0.075 0.35]); xlim([-8 28]);


%%
clear data out dayrange d x y z
end % main loop
        
figure(shallow)
orient landscape%
%print -dpdf 'Night_NEE_Surfaces_float.pdf'

figure(deep)
orient landscape
%print -dpdf 'Day_NEE_Surfaces_float.pdf'   


figure(shallow_n)
orient landscape
%print -dpdf 'Night_numer_Surfaces_float.pdf'

figure(deep_n)
orient landscape
%print -dpdf 'Day_numer_Surfaces_float.pdf'
    
    
