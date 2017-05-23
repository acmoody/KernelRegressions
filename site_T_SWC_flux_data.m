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
    function [h] =  plot( obj, varargin )

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
    h_fig1 = figure( 'Units', 'pixels', ...
                     'Position', [254 218 832 876] ); 
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
          ylabel({char(sfc_flds(fi)) ;'\theta'})
          
          if fi == 3
            xlabel('T_s_o_i_l')
          end
    end 
    suptitle(char(UNM_sites(obj.sitecode)));
  
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
         
        %ylabel(
        if fi == 1 
        title( 'all-time' );
        elseif fi == 3
            xlabel( 'T_s_o_i_l' );
        end
        
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
years = obj.clim_spaces.shallow.year_idx;

NhPanels = 3 ; 
NvPanels = 1 ; % One vertical panel for each depth

num_figs = ceil( length(years) / NhPanels ) ; 
% 
% if mod(numel(year),NhPanels)
%     num_figs = floor(9/NhPanels) + 1;
%     h1 = 4;
%     h2 = rem( 10, NhPanels );
% else
%     num_figs = numel(year)/ NhPanels;
% end % PANEL DETERMINATION
if num_figs > 1 
    if rem(length(years),NhPanels)
        n_missing_yrs = NhPanels - rem(length(years),NhPanels);
        years = min(years):1:max(years)+n_missing_yrs; 
    end
    year = reshape(years, NhPanels, num_figs);  
else
    NhPanels = length( years ); 
    year = years ;
end

    
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
        h_fig2( fi ) = figure( 'Units', 'pixels', ...
                     'Position', [129 186 1473 895] ) ;
        % Loop years
        for yi = year( 1 , fi ) : year( NhPanels , fi )  
            this_year_id = yi - year( 1 , fi ) + 1 ;
            if ~isempty(find( years == yi))
            % FOR EACH YEAR PLOT EACH DEPTH
            for di = 1:numel(depth_flds)
                print_title = false;
                refyear = min(min(year)) - 1;
                switch char(depth_flds(di))
                    case 'shallow'
                         this_cs = ...
                             obj.clim_spaces.shallow.year_clim_space( :, : , yi - refyear  );
                         print_title = true ;
                    case 'mid'
                         this_cs = ...
                             obj.clim_spaces.mid.year_clim_space( : , : , yi - refyear ) ;
                         
                    case 'deep'
                        this_cs = ...
                             obj.clim_spaces.deep.year_clim_space( : , : , yi - refyear ) ;
                end    
           % Place plot in depth di and year yi (row, column)
            this_ax = subplotrc( NvPanels, NhPanels , di,  this_year_id );
            
           
           % fprintf('v = %d h = %d di = %d yi = %d\n',v,h,di,yi)
            
            % MESH AND CONTOUR A YEAR AT A GIVEN DEPTH. ALL THE VALUES SEEM
            % TO BE THE SAME, USE ARBITRARY DEPTH
            [ T_grid, swc_grid ] = meshgrid( obj.clim_spaces.shallow.T_val, ...
                                             obj.clim_spaces.shallow.swc_val ); 
            h_temp = contourf( this_ax, ...
                      T_grid, ...
                      swc_grid, ...
                      this_cs * 100.0 ) ;
                      %this_cs.year_clim_space( :, :, this_cs_data ) * 100.0 );
            originalPosition = get(gca,'position');
            if print_title
                title( num2str( yi ) )
            end
            
            if this_year_id == 1
                  ylabel({char( depth_flds(di) );'\theta'})
            end
            h_cbar = colorbar();
            set( get( h_cbar, 'Title' ), ...
                 'String', '% days' );
            set( this_ax, 'clim', [ 0, max_frac * 100.0 ] );
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
            end % year
            end % does year exist? 
      %  fi = fi + 1;
        end  %fig

  %  suptitle( char( UNM_sites( obj.sitecode ) ) );
   
    end  % plot method
     h = [h_fig1  h_fig2(:)];
%--------------------------------------------------

end  % methods
end  % classdef
end


