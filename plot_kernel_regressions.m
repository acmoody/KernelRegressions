function h_fig = plot_kernel_regressions( kr_array, ...
                                          y_label, ...
                                          main_title, ...
                                          pdf_name, ...
                                          plot_integrated_NEE )
% PLOT_KERNEL_REGRESSIONS - 
%   

scale_term = [];

% yax_min = [ 0.025, 0.025, 0.05, 0.05, 0.05, 0.05 ];
% yax_max = [ 0.15, 0.15, 0.225, 0.225, 0.22, 0.22 ];
% xax_min = [ -8, -8, -8, -8, -8, -8 ];
% xax_max = [ 28, 28, 24, 24, 19, 19 ];   

h_fig = figure();
for this_site = 1:6
    ax = subplot( 3, 2, this_site );

    if plot_integrated_NEE
        scale_term = kr_array{ this_site }.n_count;
    end
    
    plot( kr_array{ this_site }, ...
          'ax', ax, ...
          'main_title', char( kr_array{ this_site }.sitecode ), ...
          'cbar_lab', 'NEE', ...
          'scale_term', scale_term );
    %'xlim', [ xax_min( this_site ), xax_max( this_site ) ], ...
    %      'ylim', [ yax_min( this_site ), yax_max( this_site ) ] );
end

% label the vertical axis of the middle-left plot
ax = subplot( 3, 2, 3 );
ylabel( y_label );
% label the horizontal axes of the bottom-most two plots
ax = subplot( 3, 2, 5 );
xlabel( 'Mean air temp (^oC)' );
ax = subplot( 3, 2, 6 );
xlabel( 'Mean air temp (^oC)' );

% make a 'super' title above the individual panels
suptitle( main_title );


% make an 8.5 by 11-inch PDF
set( h_fig, ...
     'Units', 'inches', ...
     'Position', [ 0, 0, 8.5, 11 ] );

fname = fullfile( getenv( 'PLOTS' ), ...
                  'KernelRegressions', ...
                  sprintf( '%s.eps', pdf_name ) );
figure_2_eps( h_fig, fname );

% for this_site = 1:1
%     % draw the figure
%     plot( kr_array{ this_site } );
%     h_fig = gcf();
%     xlabel( 'Mean air temp (^oC)' );    
%     ylabel( 'Deep SWC (cm^3 cm^-^3)' );
%    
%     title( sprintf( '%s -- %s', ...
%                     char( kr_array{ this_site }.sitecode ), ...
%                     main_title ) );
%   
%     % save the figure
%     fname = fullfile( '~/UNM/Plots/KernelRegressions', ...
%                       sprintf( '%s_deep.eps', ...
%                                char( kr_array{ this_site }.sitecode ) ) );
%     figure_2_eps( h_fig, fname );
%     %close( h_fig );
% end
