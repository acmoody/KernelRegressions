classdef site_T_SWC_flux_data
properties
    sitecode;
    surfaces;
    clim_spaces;
    note;
end

methods

%--------------------------------------------------
    function obj = site_T_SWC_flux_data( sitecode, ...
                                         sfc_shallow, ...
                                         sfc_mid, ...
                                         sfc_deep, ...
                                         clim_shallow, ...
                                         clim_mid, ...
                                         clim_deep, ...
                                         note )
    
    obj.sitecode = sitecode;
    obj.surfaces = struct( 'shallow', sfc_shallow, ...
                           'mid', sfc_mid, ...
                           'deep', sfc_deep);
    obj.clim_spaces = struct( 'shallow', clim_shallow, ...
                              'mid', clim_mid, ...
                              'deep', clim_deep);
    obj.note = note;
    end  % constructor

%--------------------------------------------------
    function [h_fig1 h_fig2] =  plot( obj, varargin )

%     h_fig1 = figure( 'Units', 'Inches', ...
%                     'Position', [ 0, 0, 24, 13 ], ...
%                     'Visible', 'on' );
    % NEW PLOTTING PLAN
    % 1 Figure for flux kernel regression and all-time climate space (for
    % each depth)
   
    % one panel for the flux kernel regression, one for each climate space year,
    % and one for the all-time climate space
%     f1h = numel( obj.clim_spaces.deep.year_idx ) + 2;
%     n_panels_v = 3; % three SWC depth panels vertically
% ------------------------------------------------------------------------
%  PLOT 1 : FLUX KERNEL REGRESSIONS AND ALL-TIME CLIMATE SPACE, ALL DEPTHS
    h_fig1 = figure
    f1h = 2 ; 
    f1v = 3 ;
    % ------------------------------------
    % PLOT THE NEE,T,SWC SURFACES 
    % loop iterators: fi - "field iterator"
    sfc_flds = fieldnames( obj.surfaces );
    for fi = 1:numel( sfc_flds )  
        this_ax = subplotrc( f1v, f1h , fi , 1 );
        plot( obj.surfaces.( sfc_flds{ fi } ), ...
              'ax', this_ax, ...
              'main_title', '', ...
              'cbar_lab', 'NEE' );
    end
    
    % -------------------------------
    % PLOT THE ALL-TIME CLIMATE SPACES
    
    all_alltime = [ obj.clim_spaces.shallow.alltime_clim_space, ...
                    obj.clim_spaces.mid.alltime_clim_space, ...
                    obj.clim_spaces.deep.alltime_clim_space ]; 
                  
    alltime_max_frac = max( reshape( all_alltime, [], 1 ) );
    
    % loop iterators: fi - "field iterator"
    clim_flds = fieldnames( obj.clim_spaces );
    for fi = 1:numel( clim_flds )
        this_cs = obj.clim_spaces.( clim_flds{ fi } );
        this_ax = subplotrc( f1v, f1h, fi, f1h );
        [ T_grid, swc_grid ] = meshgrid( this_cs.T_val, ...
                                         this_cs.swc_val ); 
        contourf( this_ax, ...
                  T_grid, ...
                  swc_grid, ...
                  this_cs.alltime_clim_space * 100.0 );
        h_cbar = colorbar();
        set( this_ax, 'clim', [ 0, alltime_max_frac * 100.0 ] );
        set( get( h_cbar, 'Title' ), ...
             'String', '% days' );
        title( 'all-time' );
        
        % draw gridlines at the SWC edges -- their spacing is contstant
        % in log space so is perhaps counterintuitive when plotted in
        % decimal space
        for i = 1:numel( this_cs.swc_val )
            h_line = refline( 0, this_cs.swc_val( i ) );
            set( h_line, 'Color', 'black', 'LineStyle', ':' );
        end
        
        fi = fi + 1;
    end % All-time climate space
    
% --------------------------------------------
% PLOT 2+, ANNUAL CLIMATE SPACES AT ALL DEPTHS

% Determine how many 3x4 figures to produce to provide better viewing of
% the yearly (seasonal?) climate spaces.
% This will remain constant over any given year.
NhPanels = 3;
year = [2007:2016];
v = numel(fieldnames(obj.clim_spaces)); % vertical panels 
if mod(numel(this_cs.year_idx),NhPanels)
    num_figs = floor(numel(this_cs.year_idx)/NhPanels) + 1;
    h1 = 4;
    h2 = rem( 10, NhPanels );
else
    num_figs = numel(this_cs.year_idx)/ NhPanels;
end % PANEL DETERMINATION
    

    
    % Normalize colors to largest fraction in climate space
    all_annual_cs = [ obj.clim_spaces.shallow.year_clim_space, ...
                      obj.clim_spaces.mid.year_clim_space, ...
                      obj.clim_spaces.deep.year_clim_space ];
    max_frac = prctile( reshape( all_annual_cs, [], 1 ) , 99.9);
    
    % plot the annual climate spaces
    % loop iterators: fi - figure iterator 
    %                 yi - year iterator"  
    %                 di - depth iterator"
    % START LOOP AT FIGURE 1
    
    depth_flds = fieldnames( obj.clim_spaces );
    for fi = 1:num_figs
        h_fig2( fi ) = figure;
        start_yr = 1+(fi-1)*4;
        end_yr = start_yr + NhPanels - 1;
            try
                year_flds = year(start_yr:end_yr);
                h = h1 ;
            catch
                year_flds = year(start_yr:end);
                h = h2 ;
            end
 %YEARS ARE SET
        for yi = 1:numel(year_flds)
            this_year_id = (year_flds) - 2007 + 1;
            % FOR EACH YEAR PLOT EACH DEPTH
            for di = 1:numel(depth_flds)
                print_title = false;
                switch char(depth_flds(di))
                    case 'shallow'
                         this_cs = ...
                             obj.clim_spaces.shallow.year_clim_space( :, : , this_year_id( yi ) );
                         print_title = true ;
                    case 'mid'
                         this_cs = ...
                             obj.clim_spaces.mid.year_clim_space( : , : , this_year_id( yi ) ) ;
                         
                    case 'deep'
                        this_cs = ...
                             obj.clim_spaces.deep.year_clim_space( : , : , this_year_id( yi ) ) ;
                end    
           % Place plot in depth di and year yi (row, column)
            this_ax = subplotrc( v, h , di,  yi ); 
            
            % MESH AND CONTOUR A YEAR AT A GIVEN DEPTH. ALL THE VALUES SEEM
            % TO BE THE SAME, USE ARBITRARY DEPTH
            [ T_grid, swc_grid ] = meshgrid( obj.clim_spaces.shallow.T_val, ...
                                             obj.clim_spaces.shallow.swc_val ); 
            h_temp = contourf( this_ax, ...
                      T_grid, ...
                      swc_grid, ...
                      this_cs * 100.0 ) ;
                      %this_cs.year_clim_space( :, :, this_cs_data ) * 100.0 );
            originalPosition = get(gca,'position')    ; 
            
            h_cbar = colorbar();
            set( get( h_cbar, 'Title' ), ...
                 'String', '% days' );
            %set( this_ax, 'clim', [ 0, max_frac * 100.0 ] );
            % RESET PLOT AXES TO ORIGINAL POSITION
            set( this_ax , 'position', originalPosition) ; 
          
%             if print_title
%             title( num2str(year_flds(this_year_id) ) );
%             end
            
            
            % draw gridlines at the SWC edges -- their spacing is contstant
            % in log space so is perhaps counterintuitive when plotted in
%             % decimal space
%             for i = 1:numel( this_cs.swc_val )
%                 h_line = refline( 0, this_cs.swc_val( i ) );
%                 set( h_line, 'Color', 'black', 'LineStyle', ':' );
%             end
        end
      %  fi = fi + 1;
    end

    suptitle( char( UNM_sites( obj.sitecode ) ) );
    
    end  % plot method
%--------------------------------------------------

end  % methods
end  % classdef
end


