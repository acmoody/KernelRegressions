function plot_T_SWC_flux_surface( ax, ...
                                  T_bin_ctrs, swc_bin_ctrs, ...
                                  flux_sfc, cmap, sitecode, bounds, ...
                                  swc_str )
% PLOT_T_SWC_FLUX_SURFACE - plots an NEE-T-SWC surface to a specified
%   axis with a specified color map


    contourf( T_bin_ctrs, exp( swc_bin_ctrs), flux_sfc ); 
    colormap( cmap );

    set(gca,'fontweight','bold','fontsize',12)
    h_cbar = colorbar;
    set( get( h_cbar, 'Title' ), 'String', 'NEE' );
    ymin = max( bounds( 1 ), min( min( exp( swc_bin_ctrs ) ) ) );
    ymax = min( bounds( 2 ), max( max( exp( swc_bin_ctrs ) ) ) );
    xmin = max( bounds( 3 ), min( min( T_bin_ctrs ) ) ); 
    xmax = min( bounds( 4 ), max( max( T_bin_ctrs ) ) );
    ylim( [ ymin, ymax ] ); 
    xlim( [ xmin, xmax ] );

    if sitecode==3
        ylabel( sprintf( '%s SWC (cm^3 cm^-^3)', swc_str ), ...
                'fontweight','bold','fontsize',14 );
    end
    if sitecode==5 || sitecode==6
        xlabel('Mean air temp (^oC)','fontweight','bold','fontsize',14);
    end