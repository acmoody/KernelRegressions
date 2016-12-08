classdef T_SWC_flux_sfc
properties
    sitecode;
    swc_val;
    T_val;
    sfc; 
    n_count;
    note;
end

methods
    
        function obj = T_SWC_flux_sfc( sitecode, swc_val, T_val, ...
                                       sfc_arg, n_count, note )
        % class constructor 
        obj.sitecode = sitecode;
        obj.swc_val = swc_val;
        obj.T_val = T_val;
        obj.sfc = sfc_arg;
        obj.n_count = n_count;
        obj.note = note;
        
        end 
    
        function plot( obj, varargin )
        % draw a contour plot of object's surface in specified axes
        
        % -----
        % define optional inputs, with defaults and typechecking
        % -----
        args = inputParser;
        args.addRequired( 'obj', @(x) isa( x, 'T_SWC_flux_sfc' ) );
        %args.addParamValue( 'ax', NaN, @isnumeric );
        args.addParamValue('ax', NaN);
        args.addParamValue( 'cmap', [], @(x) size( x, 2 ) == 3 );
        args.addParamValue( 'xlab', '', @ischar );
        args.addParamValue( 'ylab', '', @ischar );
        args.addParamValue( 'cbar_lab', '', @ischar );
        args.addParamValue( 'xlim', [], @(x) ( isnumeric( x ) & ...
                                               numel( x ) == 2 ) );
        args.addParamValue( 'ylim', [], @(x) ( isnumeric( x ) & ...
                                               numel( x ) == 2 ) );
        args.addParamValue( 'main_title', ...
                            sprintf( '%s -- %s', ...
                                     char( obj.sitecode ), ...
                                     obj.note ), ...
                            @ischar );
        args.addParamValue( 'scale_term', [], @isnumeric );
        % parse optional inputs
        args.parse( obj, varargin{ : } );
        
        % -----
        % draw the plot
        % -----
        if nargin == 1
            % if only the object itself was provided as an
            % argument, produce a basic plot with a sensible
            % set of default options.
            hfig = figure();
            ax = axes();
            cmap_YlGn = cbrewer( 'seq', 'YlGn', 9 );
            cmap_YlGn = interp1( 1:9, cmap_YlGn, linspace( 1, 9, 100 ) );
            cmap_YlGn = flipud( cmap_YlGn );
            contourf( ax, obj.T_val, obj.swc_val, obj.sfc );
            h_cbar = colorbar;
            %colormap( cmap_YlGn );
            title( args.Results.main_title );
        else
            
            %if isnan( args.Results.ax )
            if isempty(get(args.Results.ax))
                hfig = figure();
                ax = axes();
            else
                ax = args.Results.ax;
            end
            
            set( ax, 'fontweight', 'bold', 'fontsize', 12 );
            
            % if scale_term argument specified, use it to scale the surface
            if not( isempty( args.Results.scale_term ) )
                obj.sfc = obj.sfc .* args.Results.scale_term;
            end
            
            contourf( ax, obj.T_val, obj.swc_val, obj.sfc );
            
            % define the colors to use
            if not( isempty( args.Results.cmap ) )
                colormap( args.Results.cmap );
            end

            % axis limits
            if not( isempty( args.Results.xlim ) )
                set( ax, 'XLim', args.Results.xlim );
            end
            if not( isempty( args.Results.ylim ) )
                set( ax, 'YLim', args.Results.ylim );
            end
            
            % axis labels
            if not( isempty( args.Results.xlab ) )
                xlabel( args.Results.xlab );
            end
            if not( isempty( args.Results.ylab ) )
                ylabel( args.Results.ylab );
            end
            
            % draw gridlines at the SWC bin edges -- their spacing is constant
            % in log space so is perhaps counterintuitive when plotted in
            % decimal space
            swc_edges = unique( obj.swc_val );
            for i = 1:numel( swc_edges )
                h_line = refline( 0, swc_edges( i ) );
                set( h_line, 'Color', 'black', 'LineStyle', ':' );
            end
            
            % add color legend and label it
            h_cbar = colorbar();
            if not( isempty( args.Results.cbar_lab ) )
                set( get( h_cbar, 'Title' ), ...
                     'String', args.Results.cbar_lab );
            end
            
            if not( isempty( args.Results.main_title ) )
                title( args.Results.main_title );
            end
        end
        
        end   % methods
end
end