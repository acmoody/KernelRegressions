function [ andy_sfc, tim_sfc ] = scratch( sitecode )
% SCRATCH - 
%   

histogram_only = false;
kernel_regression = true;

KR_data = get_data( sitecode );
andy_sfc = andy_calc_surface( sitecode, KR_data, histogram_only );
[ tim_sfc, ~, swc_sfc, T_sfc ] = ...
    calculate_SWC_T_NEE_kernel_regression( KR_data.T, ...
                                           KR_data.swc_deep, ...
                                           KR_data.nee_day, ...
                                           histogram_only );
tim_sfc = T_SWC_flux_sfc( sitecode, exp( swc_sfc ), T_sfc, tim_sfc, ...
                          'Tim, deep, day NEE, histogram only' );

%============================================================
function KR_data = get_data( sitecode )
% GET_DATA - 

load( 'kernel_regression_parsed_data.mat' );
data = all_data{ sitecode };

data(data==-9999)=NaN;
mu2e=(1*60*30)./1000000;
mu2g=((1./1000000)*12)*60*30;

ndays=(length(data)/48);

firstday = [ 151, 90, 59, 59, 90, 120 ];
lastday = [ 272, 272, 303, 303, 303, 303 ];

for i = 1:ndays
    dayy=find(data(:,2)==i+(data(1,2)-1));
    touse=data(dayy,:);
    
    out(i,1)=nanmean(touse(:,1));       % 1 year
    out(i,2)=nanmean(touse(:,2));       % 2 day 
    out(i,3)=nanmean(touse(:,3));       % 3 hour
    out(i,4)=nanmean(touse(:,4));       % 4 decitime
    out(i,5)=nanmean(touse(:,5));       % 5 air temp
    out(i,6)=nansum(touse(:,6));           % 6 precip
    out(i,7)=nanmean(touse(:,7));       % 7 swc shallow
    out(i,8)=nanmean(touse(:,8));       % 8 swc deep
    out(i,9)=(nansum(touse(:,9))).*mu2g;   % 9 nee
    out(i,10)=(nansum(touse(:,10))).*mu2g; % 10 gpp
    out(i,11)=(nansum(touse(:,11))).*mu2g; % 11 re
    
    daytime=find(touse(:,3)>530 & touse(:,3)<1830);
    nigtime=find(touse(:,3)<600 | touse(:,3)>1800);
    
    out(i,12)=(nansum(touse(daytime,9))).*mu2g;   % 12 daytime nee
    out(i,13)=(nansum(touse(nigtime,9))).*mu2g;   % 13 nighttime nee
    out(i,14)=nanmean(touse(daytime,12));   % Mean of daytime PAR

    % remove days with precip and days following days with precip
    if i == 1
        out(i,15)=0;
    elseif out(i,6)>0
        out(i,15)=0;
    else
        out(i,15)=out(i-1,15)+1;
    end
    
end

%% Temperature response

% filter by time of year -- firstday and lastday defined above
% dayrange=find(firstday(sitecode) < out(:,4) & out(:,4)< lastday(sitecode));
dayrange=find(0 < out(:,4) & out(:,4)< 100000);

doy = out( dayrange, 2 );
swc_shallow = out( dayrange, 7 );  % shallow SWC
swc_shallow( swc_shallow < 0.001 ) = 0.001;
swc_deep = out( dayrange, 8 );  % deep SWC
swc_deep( swc_deep < 0.001 ) = 0.001;
T = out( dayrange, 5 );
nee_night = out( dayrange, 13 ); % nighttime NEE
nee_day = out( dayrange, 12 ); % nighttime NEE

swc_shallow = log( swc_shallow );
swc_deep = log( swc_deep );

KR_data = struct( 'doy', doy, ...
                  'swc_shallow', swc_shallow, ...
                  'swc_deep', swc_deep, ...
                  'T', T, ...
                  'nee_night', nee_night, ...
                  'nee_day', nee_day );


%============================================================

function andy_sfc = andy_calc_surface( sitecode, data, do_KR )
% ANDY_HISTOGRAM - 
%   

% x = out(dayrange,8); % deep SWC
% x(x<0.001)=0.001;
x = data.swc_deep;
y = data.T; %out(dayrange,5);
z = data.nee_day; %out(dayrange,12); % daytime NEE


lt=-10; % minimum temp
tr=40; % temp range
ld=log(0.07);
dr=log(0.38)-log(0.07);
tsteps=10; % number of steps
dsteps=10;
cnt=0;

swc_edges = zeros( 1, dsteps );
T_edges = zeros( 1, tsteps );

for i = 1:tsteps
    mint=(lt)+((i-1).*(tr./tsteps));
    maxt=(lt)+i.*(tr./tsteps);
    valut=(mint+maxt)./2;
    T_edges( i ) = mint;
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
        swc_edges( j ) = mind;
    end
end

[XX,YY]=meshgrid(XX2(:),YY2(:));
fluxr = flux;

if do_KR
    r=ksrmv([XX(:) YY(:)],flux(:));
    fluxr=flux;
    fluxr(:)=r.f;
end

andy_sfc = T_SWC_flux_sfc( UNM_sites( sitecode ), ...
                           exp( YY ), XX, fluxr, ...
                           'andy, deep, daytime NEE, histogram' );