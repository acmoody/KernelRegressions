%clear all
close all

load_stored_data = true;

colormap_greens = flipud( cbrewer( 'seq', 'YlGn', 100 ) );
older_KRs = cell( 6, 2 );
andy_out = cell( 6, 1 );

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

yax_min = [ 0.025, 0.025, 0.05, 0.05, 0.05, 0.05 ];
yax_max = [ 0.15, 0.15, 0.225, 0.225, 0.22, 0.22 ];
xax_min = [ -8, -8, -8, -8, -8, -8 ];
xax_max = [ 28, 28, 24, 24, 19, 19 ];

figure( 'NumberTitle', 'off', 'Name', 'shallow old' );
shallow=gcf;
figure( 'NumberTitle', 'off', 'Name', 'deep old' );
deep=gcf;
figure( 'NumberTitle', 'off', 'Name', 'deep_n old' );
deep_n=gcf;
figure( 'NumberTitle', 'off', 'Name', 'shallow_n old' );
shallow_n=gcf;

if load_stored_data
    load( 'kernel_regression_parsed_data.mat' );
else
    all_data = cell( 6, 1 );
end

for sitecode = 1:6

    % parsing takes a minutes -- option to load saved data
    if load_stored_data
        data = all_data{ sitecode };
    else
        data07 = get_kernel_regression_data( sitecode, 2007 );
        data08 = get_kernel_regression_data( sitecode, 2008 );
        data09 = get_kernel_regression_data( sitecode, 2009 );
        data10 = get_kernel_regression_data( sitecode, 2010 );
        data11 = get_kernel_regression_data( sitecode, 2011 );

        data = vertcat( data07, data08, data09, data10, data11 );
        all_data{ sitecode } = data;
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

        andy_out{ sitecode } = out;
        
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
    swc = out( dayrange, 7 );  % shallow SWC
    swc( swc < 0.001 ) = 0.001;
    T = out( dayrange, 5 );
    nee = out( dayrange, 13 ); % nighttime NEE
    
    swc = log( swc );

    % define T bins
    T_steps=10; 
    T_min = -10;
    T_max = T_min + 40;
    T_bin_edges = linspace( T_min, T_max, T_steps + 1 );
    T_bin_centers = T_bin_edges( 1:end-1 ) + ( diff( T_bin_edges ) / 2 );

    % define SWC bins
    swc_steps = 10;
    swc_min = 0.01;
    swc_max = 0.38;
    swc_bin_edges = linspace( log( swc_min ), log( swc_max ), swc_steps + 1 );
    swc_bin_centers = swc_bin_edges( 1:end-1 ) + ( diff( swc_bin_edges ) / 2 );

    % assign T, SWC observations into bins
    [ ~, T_idx ] = histc( T, T_bin_edges );
    [ ~, swc_idx ] = histc( swc, swc_bin_edges );

    % ignore observations outside of the range of T, SWC bins
    keep_idx = ( T_idx ~= 0 ) & ( swc_idx ~= 0 );
    
    % count the number of valid NEE observations in each T--SWC bin
    n_count = accumarray( [ T_idx( keep_idx ), swc_idx( keep_idx ) ], ...
                          nee( keep_idx ), ...
                          [ T_steps, swc_steps ], ...
                          @(x) sum( not( isnan( x ) ) ) );
    % calculate mean of valid NEE observations in each T--SWC bin
    flux = accumarray( [ T_idx( keep_idx ), swc_idx( keep_idx ) ], ...
                       nee( keep_idx ), ...
                       [ T_steps, swc_steps ], ...
                       @nanmean );

    [XX,YY]=meshgrid( T_bin_centers, swc_bin_centers );
    fluxr = flux;
    %remove r=ksrmv([XX(:) YY(:)],flux(:));
    %remove fluxr=flux;
    %remove fluxr(:)=r.f;
    
    %%
    left=[0.1 0.55 0.1 0.55 0.1 0.55];
    bottom=[0.7 0.7 0.4 0.4 0.1 0.1];

    figure(shallow);
    subplot('Position',[left(sitecode) bottom(sitecode) 0.4 0.25])
    
    contourf(XX,exp(YY),fluxr); 
    colormap( colormap_greens );

    if sitecode==3
        ylabel('Shallow SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
    end
    if sitecode==5 || sitecode==6
        xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
    end
    set(gca,'fontweight','bold','fontsize',12)
    h_cbar = colorbar;
    set( get( h_cbar, 'Title' ), 'String', 'NEE' );
    ymin=max(0.02,min(exp(swc)));
    ymax=min(0.31,max(exp(swc)));
    xmin=max(-8,min(T)); 
    xmax=min(28,max(T));
    ylim([ymin ymax]); 
    xlim([xmin xmax]);
    %ylim([0.02 0.31]); xlim([-8 28]);
    % ylim( [ yax_min( sitecode ), yax_max( sitecode ) ] );
    % xlim( [ xax_min( sitecode ), xax_max( sitecode ) ] );
    n_count=n_count./(sum(reshape(n_count,100,1)));

    figure(shallow_n);
    subplot('Position',[left(sitecode) bottom(sitecode) 0.4 0.25])

    contourf(XX,exp(YY),n_count.*fluxr); 
    colormap( colormap_greens );

    if sitecode==3
        ylabel('Shallow SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
    end
    if sitecode==5 || sitecode==6
        xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
    end
    set(gca,'fontweight','bold','fontsize',12)
    h_cbar = colorbar;
    set( get( h_cbar, 'Title' ), 'String', 'NEE' );

    ymin=max(0.02,min(exp(swc)));
    ymax=min(0.31,max(exp(swc)));
    xmin=max(-8,min(T)); 
    xmax=min(28,max(T));
    ylim([ymin ymax]); 
    xlim([xmin xmax]);
    %ylim([0.02 0.31]); xlim([-8 28]);
    % ylim( [ yax_min( sitecode ), yax_max( sitecode ) ] );
    % xlim( [ xax_min( sitecode ), xax_max( sitecode ) ] );

    clear x; clear y; clear z;
    %%

    x = out(dayrange,8); % deep SWC
    x(x<0.001)=0.001;
    y = out(dayrange,5);
    z = out(dayrange,12); % daytime NEE

    
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
    fluxr = flux;

    older_KRs{ sitecode, 2 } = T_SWC_flux_sfc( UNM_sites( sitecode ), ...
                                               XX, exp( YY ), fluxr, ...
                                               'deep, daytime NEE, histogram' );

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
    colormap( colormap_greens );
    
    older_KRs{ sitecode, 1 } = T_SWC_flux_sfc( UNM_sites( sitecode ), ...
                                               XX, exp( YY ), fluxr, ...
                                               'deep, daytime NEE, KR' );
    
    if sitecode==3
        ylabel('Deep SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
    end
    if sitecode==5 || sitecode==6
        xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
    end
    set(gca,'fontweight','bold','fontsize',12)
    h_cbar = colorbar;
    set( get( h_cbar, 'Title' ), 'String', 'NEE' )

    ymin=max(0.075,min(exp(x))); ymax=min(0.35,max(exp(x)));
    xmin=max(-8,min(y)); xmax=min(28,max(y));
    ylim([ymin ymax]); xlim([xmin xmax]);
    %ylim([0.075 0.35]); xlim([-8 28]);
    % ylim( [ yax_min( sitecode ), yax_max( sitecode ) ] );
    % xlim( [ xax_min( sitecode ), xax_max( sitecode ) ] );

    n_count=n_count./(sum(reshape(n_count,100,1)));
    figure(deep_n);
    subplot('Position',[left(sitecode) bottom(sitecode) 0.4 0.25])

    contourf(XX,exp(YY),n_count.*fluxr); 
    colormap( colormap_greens );
    
    if sitecode==3
        ylabel('Deep SWC (cm^3 cm^-^3)','fontweight','bold','fontsize',14)
    end
    if sitecode==5 || sitecode==6
        xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14)
    end
    set(gca,'fontweight','bold','fontsize',12)
    h_cbar = colorbar;
    set( get( h_cbar, 'Title' ), 'String', 'NEE' )

    ymin=max(0.075,min(exp(x))); ymax=min(0.35,max(exp(x)));
    xmin=max(-8,min(y)); xmax=min(28,max(y));
    ylim([ymin ymax]); xlim([xmin xmax]);
    %ylim([0.075 0.35]); xlim([-8 28]);
    % ylim( [ yax_min( sitecode ), yax_max( sitecode ) ] );
    % xlim( [ xax_min( sitecode ), xax_max( sitecode ) ] );

    %%
    clear data out dayrange d x y z
end  % site loop 

% set figure dimensions
set( shallow, 'Units', 'Inches', 'Position', [ 0, 0, 8.5, 11 ] );
set( deep, 'Units', 'Inches', 'Position', [ 0, 0, 8.5, 11 ] );
set( shallow_n, 'Units', 'Inches', 'Position', [ 0, 0, 8.5, 11 ] );
set( deep_n, 'Units', 'Inches', 'Position', [ 0, 0, 8.5, 11 ] );

% save figures to encapsulated postscript files
figure_2_eps( shallow, 'shallow.eps' );
figure_2_eps( shallow_n, 'shallow_n.eps' );
figure_2_eps( deep, 'deep.eps' );
figure_2_eps( deep_n, 'deep_n.eps' );

% if the data were parsed, save in binary format for future time savings
if not( load_stored_data )
    save( 'kernel_regression_parsed_data.mat', 'all_data' );
end
