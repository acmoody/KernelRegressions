function clim_space = calculate_climate_space( sitecode, ...
                                               T, ...
                                               SWC, ...
                                               year, ...
                                               depth )
% CALCULATE_CLIMATE_SPACE - calculate annual and across-all-data "climate
%   space".  The climate space is defined as the fraction of days within the
%   period of record within a set of defined air temperature -- soil volumetric
%   water content bins.
%
% USAGE: clim_space = calculate_climate_space( sitecode, T, SWC, year, depth )
%                                              
%
% INPUTS:
%     sitecode: UNM_sites object or integer; which site?
%     T: Nx1 numeric; air temperature
%     SWC: Nx1 numeric; soil volumetric water content
%     year: Nx1 numeric: year of each obseration
%     depth: scalar; depth of the soil water measurements (cm)
%
% OUTPUTS
%     clim_space: climate_space object
%
% (c) Timothy W. Hilton, UNM, Oct 2012
% Modified by Alex Moody, UNM, Dec 2016

[T_bin_edges, T_bin_centers, SWC_bin_edges, SWC_bin_centers] = ...
    define_kernel_regression_bins();

% place the T and SWC observations into the defined bins
[ ~, T_idx ] = histc( T, T_bin_edges );
[ ~, SWC_idx ] = histc( SWC, SWC_bin_edges );

% discard days that do not fit into any of the defined T--SWC bins
keep_idx = ( T_idx ~= 0 ) & ( SWC_idx ~= 0 );
n_keep = numel( find( keep_idx ) );

% give each day a year index for the accumarray output
[ ~, year_idx ] = ismember( year, unique( year ) );

% count days in each year within each defined T--SWC bin
annual_clim_space = ...
    accumarray( [ SWC_idx( keep_idx ), ...
                  T_idx( keep_idx ), ...
                  year_idx( keep_idx ) ], ...
                ones( n_keep, 1 ), ...
                [ numel( SWC_bin_centers ), ...
                  numel( T_bin_centers ), ...
                  numel( unique( year ) ) ] );

% count the number of observations within each year
n_days = accumarray( year_idx, ones( size( year_idx ) ) );

% if there is less than half a year of data, get rid of it!
idx = n_days < 180 ;  
annual_clim_space(:,:,idx) = [] ;
 data_years = unique(year); % get rid of that year of data, too
 data_years(idx) = [] ;
 

% divide by number of days in each year to get a fraction
for i = 1:numel( data_years )
    annual_clim_space( :, :, i ) = annual_clim_space( :, :, i ) ./ n_days( i );
end

% count days across all years within each defined T--SWC bin
all_time_clim_space = ...
    accumarray( [ SWC_idx( keep_idx ), ...
                  T_idx( keep_idx ) ], ...
                ones( n_keep, 1 ), ...
                [ numel( SWC_bin_centers ), ...
                  numel( T_bin_centers ) ] );
n_days = numel( T );  % number of days with data (all years)
all_time_clim_space = all_time_clim_space ./ n_days;

clim_space = climate_space( sitecode, ...
                            SWC_bin_centers, ...
                            T_bin_centers, ...
                            depth, ...
                            unique( year ), ...
                            annual_clim_space, ...
                            all_time_clim_space, ...
                            '' );