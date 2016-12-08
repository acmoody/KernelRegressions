clear all
%close all

load_stored_data = true;

colormap_greens = flipud( cbrewer( 'seq', 'YlGn', 100 ) );

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

figure( 'NumberTitle', 'off', 'Name', 'shallow - KR' );
shallow=gcf;
figure( 'NumberTitle', 'off', 'Name', 'deep - KR' );
deep=gcf;
figure( 'NumberTitle', 'off', 'Name', 'deep_n - KR' );
deep_n=gcf;
figure( 'NumberTitle', 'off', 'Name', 'shallow_n - KR' );
shallow_n=gcf;

figure( 'NumberTitle', 'off', 'Name', 'shallow no KR' );
shallow_extra=gcf;

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
    
    %% shallow plots
    [ flux_sfc, n_count, T_bin_ctrs, swc_bin_ctrs ] = ...
        calculate_SWC_T_NEE_kernel_regression( T, ...
                                               swc_shallow, ...
                                               nee_night, ...
                                               true);
    [ flux_sfc_nokr, n_count, T_bin_ctrs, swc_bin_ctrs ] = ...
        calculate_SWC_T_NEE_kernel_regression( T, ...
                                               swc_shallow, ...
                                               nee_night, ...
                                               false );
    shallow_plot_bounds = [ 0.02, 0.31, -8, 28 ];
    deep_plot_bounds = [ 0.075, 0.35, -8, 28 ];
    figure(shallow);
    ax = subplot( 3, 2, sitecode );
    plot_T_SWC_flux_surface( ax, ...
                             T_bin_ctrs, swc_bin_ctrs, ...
                             flux_sfc, colormap_greens, sitecode, ...
                             shallow_plot_bounds, ...
                             'Shallow' );
    
    figure(shallow_extra);
    ax = subplot( 3, 2, sitecode );
    plot_T_SWC_flux_surface( ax, ...
                             T_bin_ctrs, swc_bin_ctrs, ...
                             flux_sfc_nokr, colormap_greens, sitecode, ...
                             shallow_plot_bounds, ...
                             'Shallow' );

    figure(shallow_n);
    ax = subplot( 3, 2, sitecode );
    n_count = n_count ./ ( sum( reshape( n_count, [], 1 ) ) );
    flux_sfc = flux_sfc .* n_count;
    plot_T_SWC_flux_surface( ax, ...
                             T_bin_ctrs, swc_bin_ctrs, ...
                             flux_sfc, colormap_greens, sitecode, ...
                             shallow_plot_bounds, ...
                             'Shallow' );

    %% deep plots
    [ flux_sfc, n_count, T_bin_ctrs, swc_bin_ctrs ] = ...
        calculate_SWC_T_NEE_kernel_regression( T, ...
                                               swc_deep, ...
                                               nee_day, ...
                                               true);
    
    figure(deep);
    ax = subplot( 3, 2, sitecode );
    plot_T_SWC_flux_surface( ax, ...
                             T_bin_ctrs, swc_bin_ctrs, ...
                             flux_sfc, colormap_greens, sitecode, ...
                             deep_plot_bounds, ...
                             'Deep' );

    figure(deep_n);
    ax = subplot( 3, 2, sitecode );
    n_count = n_count ./ ( sum( reshape( n_count, [], 1 ) ) );
    flux_sfc = flux_sfc .* n_count;
    plot_T_SWC_flux_surface( ax, ...
                             T_bin_ctrs, swc_bin_ctrs, ...
                             flux_sfc, colormap_greens, sitecode, ...
                             deep_plot_bounds, ...
                             'Deep' );

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
