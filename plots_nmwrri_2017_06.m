% ET PLOTS FOR NMWRRI JUNE 2017 ET CONFERENCE
% The kernel regression plotting is a mess. Just laod surface data and do
% plotting here. 

load('surface_data.mat')
load('C:\Code\KernelRegressions\kernel_regression_parsed_data.mat')
% Surface_data.mat is a terribly complicated cell. 

% --------- Plot ET surface for 6 Ameriflux core sites 


h_fig = figure( 'Units', 'Inches', ...
                'Position', [ 0, 0, 24, 13 ], ...
                'Visible', 'on' );

% one panel for the flux kernel regression for each site
n_panels_h = 3;
n_panels_v = 2; % four SWC depth panels vertically

%pal = colormap( cbrewer( 'seq', 'Blues', 5 ) );
%ET_cmap = [ interp1( 1:5, pal, linspace( 1, 5, 100 ) ) ];
pal = colormap('jet');
pal = flipud(pal);


% plot the ET--T--SWC surfaces

for i = 1:6
    h{i} = subplot( n_panels_v, n_panels_h, i );
    plot( surface_data{i,1}.surfaces.shallow, ...
          'ax', h{i}, ...
          'cmap',pal,...
          'cbar_lab', 'ET mm' );
end

set([h{:}],'clim',[0.5 3.5]);

%%
% ----------------------- Plot WY Bar plots ----------------
% Add hydroyear and hydro DOY first
for i = 1:numel(all_data)
    ts = all_data{i,1}.TIMESTAMP;
    [all_data{i,1}.year, all_data{i,1}.month , all_data{i,1}.day,~,~,~] = ...
        datevec(ts);
    yearvec = unique(all_data{i,1}.year);
    hydroyear = zeros(height(all_data{i,1}),1);
    season = zeros(height(all_data{i,1}),1);
    for j = 1:length(yearvec)
        wy_idx = find( ts <= datenum( yearvec(j) , 9 , 30) & ...
                       ts >= datenum( yearvec(j) - 1 , 10, 1) );
       hydroyear(wy_idx) = yearvec(j) - 1  ;
    end
    all_data{i,1}.hydroyear = hydroyear;
    all_data{i,1}.hDOY =  all_data{i,1}.TIMESTAMP - datenum(  all_data{i,1}.hydroyear, 10, 1 ) +1;
    cold_idx = find(all_data{i,1}.hDOY <= 182 );
    season(cold_idx,:) = 1; % COLD SEASON
    spring_idx = find(all_data{i,1}.hDOY > 182 & all_data{i,1}.hDOY <= 273 );
    season(spring_idx) = 2; % SPRING
    monsoon_idx = find(all_data{i,1}.hDOY > 273 & all_data{i,1}.hDOY <= 366 );
    season(monsoon_idx) = 3; % MONSOON
    all_data{i,1}.season = season;
    
    % Cumulative Fluxes
    [uniqueWY,idxToUnique,idxFromUniqueBackToAll] = unique(all_data{i,1}.hydroyear);
    % Accumulatve Precip and ET values over water year
    cumulativeET = accumarray(idxFromUniqueBackToAll,all_data{i,1}.ET_mm_dayint,[],@(x) {cumsum(x,'omitnan')});
    all_data{i,1}.cumET = vertcat(cumulativeET{:});
    cumulativePPT = accumarray(idxFromUniqueBackToAll,all_data{i,1}.PRECIP,[],@(x) {cumsum(x,'omitnan')});
    all_data{i,1}.cumPPT = vertcat(cumulativePPT{:});
    
end

%  Plot Yearly Cumulative values
%%

for j = 1:length(uniqueWY)
    h_fig = figure;
   % idx = find(all_data{i,1}.hydroyear == uniqueWY(j) );
    for i = 1:numel(all_data)
        idx = find(all_data{i,1}.hydroyear == uniqueWY(j) );
        this_ts = datenum(all_data{i,1}.hydroyear(idx),1,0)+all_data{i,1}.hDOY(idx);
        this_ts = all_data{i,1}.TIMESTAMP(idx); 
        plot( this_ts,[all_data{i, 1}.cumPPT(idx) - all_data{i,1}.cumET(idx)] ,'LineWidth',2);
        hold on
    end
    hold off
    legend('GLand','SLand','JSav','PJ','PPine','MCon','Location','Best')
    datetick('x','mmm-yy')
    ylabel('P - ET [mm]')
    xlabel('Date')
    title(sprintf('Water Year %d',uniqueWY(j)))
    saveas(h_fig,...
        fullfile('C:\Research_Flux_Towers\Plots\NMWRRI_ET',sprintf('precip_less_et_%d.png',uniqueWY(j))))
end
% ------------ Mean Annual ET/P for each year at each site
%%

plot_mean_et_with_error


