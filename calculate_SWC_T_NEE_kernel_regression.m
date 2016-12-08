function [ flux_sfc, n_count, swc_grid, T_grid ] = ...
    calculate_SWC_T_NEE_kernel_regression( T, swc, flux, do_kernel_reg )
% CALCULATE_SWC_T_NEE_KERNEL_REGRESSION - uses kernel regression to
%   fit a a surface to flux as a function of soil water content and
%   temperature.  T and SWC are grouped into ten bins and flux
%   tabulated within each T--SWC bin.  A 10 x 10 surface is then
%   fit to the tabulated flux.
%
% USAGE
%     [ flux_sfc, n_count, XX, YY ] = ...
%         calculate_SWC_T_NEE_kernel_regression( T, swc, flux )
%
% INPUTS:
%     T: N-element numeric vector; soil temperature
%     swc: N-element numeric vector; soil water content
%     flux:  N-element numeric vector; carbon flux (e.g. GPP)
%
% OUTPUTS
%     flux_sfc: the 10 x 10 surface fit to the flux
%     n_count: 10 x 10 matrix counting the number of flux
%         observations in each T--SWC bin.
%     swc_grid: the 10x10 grid of SWCs defining the surface 
%     T_grid: the 10x10 grid of temperatures defining the surface 
%
% (c) Timothy W. Hilton, UNM, Sep 2012

[T_bin_edges, T_bin_centers, swc_bin_edges, swc_bin_centers] = ...
    define_kernel_regression_bins();

% assign T, SWC observations into bins
[ ~, T_idx ] = histc( T, T_bin_edges );
[ ~, swc_idx ] = histc( swc, swc_bin_edges );

% ignore observations outside of the range of T, SWC bins
keep_idx = ( T_idx ~= 0 ) & ( swc_idx ~= 0 );

% count the number of valid NEE observations in each T--SWC bin
fillval = 0;
n_count = accumarray( [ swc_idx( keep_idx ), T_idx( keep_idx ) ], ...
                      flux( keep_idx ), ...
                      [ numel( swc_bin_centers ), numel( T_bin_centers ) ], ...
                      @(x) sum( not( isnan( x ) ) ), ...
                      fillval );
% calculate mean of valid NEE observations in each T--SWC bin
fillval = NaN;
mean_flux = accumarray( [ swc_idx( keep_idx ), T_idx( keep_idx ) ], ...
                        flux( keep_idx ), ...
                        [ numel( swc_bin_centers ), ...
                          numel( T_bin_centers ) ], ...
                        @nanmean, ...
                        fillval );

[ T_grid, swc_grid ] = meshgrid( T_bin_centers, swc_bin_centers );
if do_kernel_reg
    r = ksrmv( [ swc_grid(:), T_grid(:) ], mean_flux(:) );
    flux_sfc = reshape( r.f, 10, 10 );
else
    fprintf( [ 'NOT PERFORMING KERNEL REGRESSION - ', ...
               'returning plain histogram\n' ] );
    flux_sfc = mean_flux;
end

% transform SWC grid back out of log space before returning
% keyboard
% swc_bin_centers = exp( swc_bin_centers );
% swc_grid = exp( swc_grid );
