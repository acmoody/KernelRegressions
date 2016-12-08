function [T_bin_edges, T_bin_centers, swc_bin_edges, swc_bin_centers] = ...
    define_kernel_regression_bins()
% DEFINE_KERNEL_REGRESSION_BINS - defines the bins for grouping soil water
%   content and air temperature observations to use in SWC--T--flux surface
%   estimations (i.e. kernel regressions).  
%
% This is defined in a function so that the bins can be defined in one place,
% but accessed in multiple places (e.g. kernel regression calculation,
% climate space calculation)

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
swc_bin_edges = linspace( 0, 0.40, 11 );
% swc_bin_edges = [ -2.659260, -2.490092, -2.320925, -2.151757, -1.982590, ...
%                     -1.813422, -1.644254, -1.475087, -1.305919, -1.136752, ...
%                     -0.967584 ];
swc_bin_centers = swc_bin_edges( 1:end-1 ) + ( diff( swc_bin_edges ) / 2 );

%swc_bin_edges = exp( swc_bin_edges );
%swc_bin_centers = exp( swc_bin_centers );